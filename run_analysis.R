

METHODS <- c("full", "fixed_baseline", "independent_baseline", "ignore_outcome_dependence")
TARGETS <- c("early_death", "early_death_or_deterioration", "logRR_contrast")
CONTRASTS <- rbind(early_death = c(1,0), early_death_or_deterioration = c(0,1), logRR_contrast = c(1,-1))
Z95 <- qnorm(.975)
sym <- function(x) (x+t(x))/2


check_covariance <- function(S, positive = FALSE, context = "covariance") {
  if (!is.matrix(S) || nrow(S) != ncol(S) || any(!is.finite(S)))
    stop(context, ": invalid dimensions or entries.")
  scale <- max(1, max(abs(S)))
  if (max(abs(S - t(S))) > 1e-10 * scale) stop(context, ": not symmetric.")
  eig <- eigen(sym(S), symmetric = TRUE, only.values = TRUE)$values
  if (min(eig) < -1e-10 * scale || (positive && min(eig) <= 0))
    stop(context, ": covariance is not ", if (positive) "positive definite." else "positive semidefinite.")
  invisible(TRUE)
}

binary_bounds <- function(p1, p2) {
  if (any(!is.finite(c(p1, p2))) || any(c(p1, p2) <= 0 | c(p1, p2) >= 1))
    stop("Correlation bounds require marginal probabilities strictly in (0,1).")
  scale <- sqrt(p1 * (1 - p1) * p2 * (1 - p2))
  c(lower = (max(0, p1 + p2 - 1) - p1 * p2) / scale,
    upper = (min(p1, p2) - p1 * p2) / scale)
}

# Order p = (D_t, DD_t, D_c, DD_c), arms independent. q is P(D and DD).
# Under verified same-population nesting: joint_t = p[1], joint_c = p[3].

risk_covariance <- function(p, nt, nc, joint_t = p[1], joint_c = p[3]) {
  if (length(p) != 4L || any(!is.finite(p)) || any(p <= 0 | p >= 1) ||
      length(nt) != 1L || length(nc) != 1L ||
      any(!is.finite(c(nt, nc))) || min(nt, nc) <= 0)
    stop("Invalid risks or sample sizes.")
  C <- matrix(0, 4, 4)
  for (a in 1:2) {
    ix <- if (a == 1L) 1:2 else 3:4
    n <- if (a == 1L) nt else nc
    q <- if (a == 1L) joint_t else joint_c
    v <- p[ix]
    if (!is.finite(q) || q < max(0, sum(v) - 1) - 1e-12 ||
        q > min(v) + 1e-12) stop("Joint probability violates binary marginal bounds.")
    block <- diag(v * (1 - v) / n, 2)
    block[1, 2] <- block[2, 1] <- (q - prod(v)) / n
    C[ix, ix] <- block
  }
  check_covariance(C)
  C
}

# U = (reported effect 1, reported effect 2, baseline risk 1, baseline risk 2).
# convert[j] = TRUE: reported log OR is converted; FALSE: log RR is already reported.
# All covariance blocks are constructed from the SAME underlying risk vector.

transform_effects <- function(p, C, convert = c(TRUE, TRUE)) {
  if (length(p) != 4L || any(!is.finite(p)) || any(p <= 0 | p >= 1) ||
      !all(dim(C) == c(4L, 4L)) || length(convert) != 2L || anyNA(convert))
    stop("Invalid transformation inputs.")
  check_covariance(C)
  pt <- p[1:2]; pc <- p[3:4]
  L <- stats::qlogis(pt) - stats::qlogis(pc)
  direct <- log(pt) - log(pc)
  D <- cbind(diag(1 / pt, 2), -diag(1 / pc, 2))
  H <- matrix(0, 4, 4)
  H[1:2, 1:2] <- diag(1 / (pt * (1 - pt)), 2)
  H[1:2, 3:4] <- -diag(1 / (pc * (1 - pc)), 2)
  H[3:4, 3:4] <- diag(2)
  J <- cbind(diag(1 - pt, 2), diag((pc - pt) / (pc * (1 - pc)), 2))
  y <- L - log1p(pc * expm1(L))
  for (j in 1:2) if (!convert[j]) {
    H[j, ] <- D[j, ]
    J[j, ] <- 0; J[j, j] <- 1
    y[j] <- direct[j]
  }
  omega <- sym(H %*% C %*% t(H))
  full <- sym(J %*% omega %*% t(J))
  JE <- J[, 1:2, drop = FALSE]; JB <- J[, 3:4, drop = FALSE]
  fixed <- sym(JE %*% omega[1:2, 1:2] %*% t(JE))
  independent <- sym(fixed + JB %*% omega[3:4, 3:4] %*% t(JB))
  list(y = y, L = L, direct_y = direct, direct_cov = sym(D %*% C %*% t(D)),
       omega = omega, H = H, J = J, D = D,
       S = list(full = full, fixed_baseline = fixed,
                independent_baseline = independent,
                ignore_outcome_dependence = diag(diag(full), 2)))
}

summarize_fit <- function(mu, V) {
  v <- diag(CONTRASTS %*% V %*% t(CONTRASTS))
  if (any(!is.finite(v)) || any(v <= 0)) stop("Nonpositive summary variance.")
  est <- drop(CONTRASTS %*% mu); se <- sqrt(v)
  data.frame(target = TARGETS, estimate = est, se = se,
             lower = est - Z95 * se, upper = est + Z95 * se, row.names = NULL)
}

gls_fit <- function(y, S) {
  if (!is.matrix(y) || ncol(y) != 2L || length(S) != nrow(y) || any(!is.finite(y)))
    stop("GLS input dimensions or values are invalid.")
  W <- lapply(seq_along(S), function(i) {
    check_covariance(S[[i]], positive = TRUE, context = paste("study", i))
    solve(S[[i]])
  })
  V <- sym(solve(Reduce(`+`, W)))
  rhs <- Reduce(`+`, lapply(seq_len(nrow(y)), function(i) W[[i]] %*% y[i, ]))
  mu <- drop(V %*% rhs)
  list(mu = mu, V = V, summary = summarize_fit(mu, V))
}

# REML objective (up to a constant): 1/2 [sum log|S_i+Psi| +
# log|sum (S_i+Psi)^(-1)| + sum (y_i-mu)'(S_i+Psi)^(-1)(y_i-mu)].

nested_joint_counts <- function(n_D, n_DD, d, dd, allow_fractional = FALSE) {
  x <- c(n_D, n_DD, d, dd)
  if (length(x) != 4L || any(!is.finite(x)) || n_D <= 0 || n_DD <= 0)
    stop("Nested reconstruction needs four finite scalars and positive totals.")
  if (!allow_fractional && any(abs(x - round(x)) > 1e-10))
    stop("Observed nested counts must be integers.")
  if (abs(n_D - n_DD) > 1e-10)
    stop("Unequal denominators: do not subtract marginal counts to invent a joint table.")
  if (d < 0 || dd < d || dd > n_D)
    stop("Counts must satisfy 0 <= D <= DD <= n.")
  c(D_and_DD = d, D_without_DD = 0, DD_without_D = dd - d, neither = n_D - dd)
}

nested_arm <- function(n, d, dd) {
  cells <- nested_joint_counts(n, n, d, dd)
  p <- c(D = d / n, DD = dd / n)
  C <- matrix(c(p[1] * (1 - p[1]), p[1] * (1 - p[2]),
                p[1] * (1 - p[2]), p[2] * (1 - p[2])), 2, 2) / n
  dimnames(C) <- list(c("D", "DD"), c("D", "DD"))
  q <- unname(cells[c("D_and_DD", "DD_without_D", "neither")]) / n
  M <- rbind(c(1, 0, 0), c(1, 1, 0))
  C_mult <- M %*% (diag(q, 3) - tcrossprod(q)) %*% t(M) / n
  rho <- if (all(p > 0 & p < 1))
    as.numeric(sqrt(p[1] * (1 - p[2]) / (p[2] * (1 - p[1])))) else NA_real_
  list(cells = cells, probabilities = p, covariance = C,
       multinomial_covariance = C_mult, rho = rho)
}

nested_risk_covariance <- function(p, nt, nc) {
  if (length(p) != 4L || any(!is.finite(p)) ||
      p[1] > p[2] || p[3] > p[4]) stop("Invalid nested risk ordering.")
  risk_covariance(p, nt, nc, joint_t = p[1], joint_c = p[3])
}

nested_effects <- function(nt, dt, ddt, nc, dc, ddc) {
  T <- nested_arm(nt, dt, ddt); C0 <- nested_arm(nc, dc, ddc)
  p <- c(unname(T$probabilities), unname(C0$probabilities))
  if (any(p <= 0 | p >= 1))
    stop("A marginal event/non-event count is zero: log effects need a declared handling rule.")
  C <- nested_risk_covariance(p, nt, nc)
  rr <- transform_effects(p, C)
  G_OR <- cbind(diag(1 / (p[1:2] * (1 - p[1:2])), 2),
                -diag(1 / (p[3:4] * (1 - p[3:4])), 2))
  G_RD <- cbind(diag(2), -diag(2))
  S_OR <- sym(G_OR %*% C %*% t(G_OR))
  S_RR <- rr$direct_cov
  S_RD <- sym(G_RD %*% C %*% t(G_RD))
  cov_OR_closed <- 1 / (nt * p[2] * (1 - p[1])) +
                   1 / (nc * p[4] * (1 - p[3]))
  cov_RR_closed <- (1 - p[2]) / (nt * p[2]) + (1 - p[4]) / (nc * p[4])
  cov_RD_closed <- p[1] * (1 - p[2]) / nt + p[3] * (1 - p[4]) / nc
  wh10 <- function(n, d, dd) (1 / dd + 1 / (n - dd)) * (1 - (dd - d) / (n - d))
  cov_OR_wh10 <- wh10(nt, dt, ddt) + wh10(nc, dc, ddc)
  A <- rbind(c(1, -1), c(0, 1))
  decorrelated <- sym(A %*% S_RR %*% t(A))
  checks <- c(multinomial = max(abs(T$covariance - T$multinomial_covariance),
                                 abs(C0$covariance - C0$multinomial_covariance)),
              logOR_closed = abs(S_OR[1, 2] - cov_OR_closed),
              logOR_WeiHiggins_complete = abs(S_OR[1, 2] - cov_OR_wh10),
              logRR_closed = abs(S_RR[1, 2] - cov_RR_closed),
              logRR_cov_equals_var_DD = abs(S_RR[1, 2] - S_RR[2, 2]),
              logRR_contrast_variance = abs(sum(S_RR * outer(c(1, -1), c(1, -1))) -
                                           (S_RR[1, 1] - S_RR[2, 2])),
              logRR_decorrelation = abs(decorrelated[1, 2]),
              RD_closed = abs(S_RD[1, 2] - cov_RD_closed),
              transformed_covariance = max(abs(rr$S$full - S_RR)))
  list(treatment = T, control = C0, risks = p, y_RR = rr$direct_y,
       y_OR = rr$L, y_RD = p[1:2] - p[3:4], S_RR = S_RR,
       S_OR = S_OR, S_RD = S_RD, transformed = rr,
       contrast_DD_covariance = decorrelated, checks = checks,
       singular_logRR = min(eigen(S_RR, symmetric = TRUE, only.values = TRUE)$values) <= 1e-12)
}

correct_joint_counts <- function(ct, cc) {
  if (length(ct) != 3L || length(cc) != 3L) stop("Expected three joint categories per arm.")
  dims <- dim(ct[[1]])
  if (length(dims) != 2L) stop("Counts must be B x K matrices.")
  for (m in c(ct, cc)) if (!identical(dim(m), dims) || any(!is.finite(m)) || any(m < 0))
    stop("Invalid joint count array.")
  zero <- Reduce(`|`, lapply(c(ct, cc), function(x) x == 0))
  list(ct = lapply(ct, function(x) x + .5 * zero),
       cc = lapply(cc, function(x) x + .5 * zero), adjusted = zero)
}

batch_quantities <- function(ct, cc, reporting) {
  nt <- Reduce(`+`, ct); nc <- Reduce(`+`, cc)
  t1 <- ct[[1]] / nt; t2 <- (ct[[1]] + ct[[2]]) / nt
  b1 <- cc[[1]] / nc; b2 <- (cc[[1]] + cc[[2]]) / nc
  B <- nrow(t1); K <- ncol(t1)
  pattern <- if (reporting == "all_OR") matrix(TRUE, K, 2) else
    outer(seq_len(K) - 1L, 0:1, function(i, j) (i + j) %% 2L == 0L)
  conv1 <- matrix(rep(pattern[, 1], each = B), B, K)
  conv2 <- matrix(rep(pattern[, 2], each = B), B, K)
  vt11 <- t1 * (1 - t1) / nt; vt22 <- t2 * (1 - t2) / nt
  vc11 <- b1 * (1 - b1) / nc; vc22 <- b2 * (1 - b2) / nc
  vt12 <- t1 * (1 - t2) / nt; vc12 <- b1 * (1 - b2) / nc
  y1 <- log(t1) - log(b1); y2 <- log(t2) - log(b2)
  make <- function(c1, c2) list(s11 = vt11 / t1^2 + vc11 * c1^2,
                               s22 = vt22 / t2^2 + vc22 * c2^2,
                               s12 = vt12 / (t1 * t2) + vc12 * c1 * c2)
  full <- make(1 / b1, 1 / b2)
  a1 <- ifelse(conv1, (1 - t1) / (b1 * (1 - b1)), 1 / b1)
  a2 <- ifelse(conv2, (1 - t2) / (b2 * (1 - b2)), 1 / b2)
  fixed <- make(a1, a2)
  h1 <- ifelse(conv1, (b1 - t1) / (b1 * (1 - b1)), 0)
  h2 <- ifelse(conv2, (b2 - t2) / (b2 * (1 - b2)), 0)
  independent <- list(s11 = fixed$s11 + vc11 * h1^2,
                      s22 = fixed$s22 + vc22 * h2^2,
                      s12 = fixed$s12 + vc12 * h1 * h2)
  ignore <- full; ignore$s12 <- matrix(0, B, K)
  list(y1 = y1, y2 = y2, pattern = pattern, t1 = t1, t2 = t2,
       b1 = b1, b2 = b2, nt = nt, nc = nc,
       S = list(full = full, fixed_baseline = fixed,
                independent_baseline = independent, ignore_outcome_dependence = ignore))
}

batch_gls <- function(y1, y2, S) {
  dets <- S$s11 * S$s22 - S$s12^2
  okstudy <- is.finite(dets) & dets > 0 & S$s11 > 0 & S$s22 > 0 &
             is.finite(y1) & is.finite(y2)
  valid <- rowSums(!okstudy) == 0L
  dd <- dets; dd[!okstudy] <- NA_real_
  w11 <- S$s22 / dd; w22 <- S$s11 / dd; w12 <- -S$s12 / dd
  a11 <- rowSums(w11); a22 <- rowSums(w22); a12 <- rowSums(w12)
  detA <- a11 * a22 - a12^2
  valid <- valid & is.finite(detA) & detA > 0
  v11 <- a22 / detA; v22 <- a11 / detA; v12 <- -a12 / detA
  r1 <- rowSums(w11 * y1 + w12 * y2)
  r2 <- rowSums(w12 * y1 + w22 * y2)
  mu1 <- v11 * r1 + v12 * r2; mu2 <- v12 * r1 + v22 * r2
  vars <- cbind(v11, v22, v11 + v22 - 2 * v12)
  est <- cbind(mu1, mu2, mu1 - mu2)
  valid <- valid & rowSums(!is.finite(vars) | vars <= 0 | !is.finite(est)) == 0L
  est[!valid, ] <- NA_real_; vars[!valid, ] <- NA_real_
  list(estimate = est, se = sqrt(vars), valid = valid)
}

validate_batch <- function(q, ct, cc) {
  maxerr <- 0; K <- ncol(q$y1)
  for (j in seq_len(K)) {
    p <- c(q$t1[1, j], q$t2[1, j], q$b1[1, j], q$b2[1, j])
    o <- transform_effects(p, risk_covariance(p, q$nt[1, j], q$nc[1, j]), q$pattern[j, ])
    for (m in METHODS) {
      s <- q$S[[m]]
      mat <- matrix(c(s$s11[1, j], s$s12[1, j], s$s12[1, j], s$s22[1, j]), 2, 2)
      maxerr <- max(maxerr, abs(o$S[[m]] - mat))
    }
  }
  ss <- lapply(seq_len(K), function(j) {
    s <- q$S$full
    matrix(c(s$s11[1, j], s$s12[1, j], s$s12[1, j], s$s22[1, j]), 2, 2)
  })
  f <- gls_fit(cbind(q$y1[1, ], q$y2[1, ]), ss)$summary
  b <- batch_gls(q$y1, q$y2, q$S$full)
  maxerr <- max(maxerr, abs(f$estimate - b$estimate[1, ]), abs(f$se - b$se[1, ]))
  if (maxerr > 1e-9) stop("Batch/scalar validation failed.")
  maxerr
}

simulation_scenarios <- function() data.frame(
  scenario = 1:4, regime = c("moderate", "moderate", "sparse", "sparse"),
  n_per_arm = c(100L, 500L, 100L, 500L),
  p_D = c(.10, .10, .02, .02), p_DD = c(.30, .30, .10, .10))

load_replay <- function(root, B) {
  path <- file.path(root, "inputs", "simulation_counts.csv.gz")
  if (!file.exists(path)) stop("Replay counts are missing. Use the full ZIP or --simulation generate.")
  checksum_path <- file.path(root, "inputs", "simulation_counts.md5")
  if (file.exists(checksum_path)) {
    expected <- trimws(readLines(checksum_path, warn = FALSE)[1])
    if (unname(tools::md5sum(path)) != expected) stop("Replay input checksum mismatch.")
  }
  f <- utils::read.csv(path, stringsAsFactors = FALSE)
  required <- c("scenario", "replicate", "study", "c_D", "c_DD_only", "c_neither",
                "t_D", "t_DD_only", "t_neither")
  if (!all(required %in% names(f))) stop("Replay file has unexpected columns.")
  if (any(!is.finite(as.matrix(f[, required]))) || any(as.matrix(f[, required]) < 0))
    stop("Replay counts or identifiers are invalid.")
  if (B > min(tapply(f$replicate, f$scenario, max)))
    stop("Replay has 2,000 repetitions per scenario. Use --simulation generate for a larger run.")
  f <- f[f$replicate <= B, ]
  f <- f[order(f$scenario, f$replicate, f$study), ]
  if (anyDuplicated(f[, c("scenario", "replicate", "study")])) stop("Duplicate replay rows.")
  f
}

scenario_counts <- function(cfg, B, K, replay, seed) {
  pc <- c(cfg$p_D, cfg$p_DD); pt <- pc * c(.8, .9)
  category <- function(p) c(p[1], p[2] - p[1], 1 - p[2])
  if (any(category(pt) <= 0) || any(category(pc) <= 0)) stop("Invalid nested simulation probabilities.")
  if (is.null(replay)) {
    set.seed(as.integer((as.double(seed) + cfg$scenario) %% .Machine$integer.max))
    cmat <- t(stats::rmultinom(B * K, cfg$n_per_arm, category(pc)))
    tmat <- t(stats::rmultinom(B * K, cfg$n_per_arm, category(pt)))
  } else {
    x <- replay[replay$scenario == cfg$scenario, ]
    if (nrow(x) != B * K || !identical(as.integer(x$replicate), rep(seq_len(B), each = K)) ||
        !identical(as.integer(x$study), rep(seq_len(K), times = B))) stop("Incomplete replay scenario.")
    cmat <- as.matrix(x[, c("c_D", "c_DD_only", "c_neither")])
    tmat <- as.matrix(x[, c("t_D", "t_DD_only", "t_neither")])
  }
  if (any(rowSums(cmat) != cfg$n_per_arm) || any(rowSums(tmat) != cfg$n_per_arm))
    stop("Joint counts do not sum to the planned arm size.")
  list(cc = lapply(1:3, function(j) matrix(cmat[, j], B, K, byrow = TRUE)),
       ct = lapply(1:3, function(j) matrix(tmat[, j], B, K, byrow = TRUE)))
}

simulation_metrics <- function(fit, reference, truth, cfg, reporting, method, B) {
  rows <- vector("list", 3L)
  for (j in 1:3) {
    ok <- fit$valid & is.finite(fit$estimate[, j]) & is.finite(fit$se[, j])
    n <- sum(ok); est <- fit$estimate[ok, j]; se <- fit$se[ok, j]
    err <- est - truth[j]
    covered <- rep(NA, B)
    covered[ok] <- abs(err) <= Z95 * se
    ref_ok <- reference$valid & is.finite(reference$estimate[, j]) & is.finite(reference$se[, j])
    paired <- ok & ref_ok
    delta <- as.numeric(covered[paired]) - as.numeric(
      abs(reference$estimate[paired, j] - truth[j]) <= Z95 * reference$se[paired, j])
    coverage <- if (n) mean(covered[ok]) else NA_real_
    rmse <- if (n) sqrt(mean(err^2)) else NA_real_
    emp <- if (n > 1L) stats::sd(est) else NA_real_
    rows[[j]] <- data.frame(scenario = cfg$scenario, regime = cfg$regime,
      n_per_arm = cfg$n_per_arm, reporting = reporting, method = method,
      target = TARGETS[j], reps = B, successful_reps = n, truth = truth[j],
      mean_estimate = if (n) mean(est) else NA_real_, bias = if (n) mean(err) else NA_real_,
      bias_mcse = emp / sqrt(n), rmse = rmse,
      rmse_mcse = if (n > 1L && rmse > 0) stats::sd(err^2) / (2 * rmse * sqrt(n)) else NA_real_,
      empirical_se = emp, mean_reported_se = if (n) mean(se) else NA_real_,
      coverage = coverage, coverage_mcse = if (n) sqrt(coverage * (1 - coverage) / n) else NA_real_,
      unconditional_coverage = sum(covered[ok]) / B,
      mean_interval_width = if (n) 2 * Z95 * mean(se) else NA_real_,
      width_mcse = if (n > 1L) 2 * Z95 * stats::sd(se) / sqrt(n) else NA_real_,
      coverage_difference_vs_full = if (length(delta)) mean(delta) else NA_real_,
      paired_difference_mcse = if (length(delta) > 1L) stats::sd(delta) / sqrt(length(delta)) else NA_real_,
      paired_successful_reps = length(delta), failure_rate = 1 - n / B)
  }
  do.call(rbind, rows)
}



make_audit <- function(d) {
  eq <- d$nt_D == d$nt_DD & d$nc_D == d$nc_DD
  ordered <- d$st_D <= d$st_DD & d$sc_D <= d$sc_DD
  valid <- d$nt_D > 0 & d$nc_D > 0 & d$nt_DD > 0 & d$nc_DD > 0 &
    d$st_D >= 0 & d$sc_D >= 0 & d$st_DD >= 0 & d$sc_DD >= 0 &
    d$st_D <= d$nt_D & d$sc_D <= d$nc_D &
    d$st_DD <= d$nt_DD & d$sc_DD <= d$nc_DD
  is_barer <- grepl("^Barer", d$studyID)
  data.frame(
    source_row = d$source_row, study = d$studyID,
    n_t_early_death = d$nt_D, n_c_early_death = d$nc_D,
    n_t_composite = d$nt_DD, n_c_composite = d$nc_DD,
    events_t_death = d$st_D, events_c_death = d$sc_D,
    events_t_composite = d$st_DD, events_c_composite = d$sc_DD,
    equal_denominators = eq, event_counts_ordered = ordered,
    counts_valid = valid,
    numeric_candidate = eq & ordered & valid,
    barer_shared_control = is_barer,
    included_in_primary_conditional_analysis = eq & ordered & valid & !is_barer,
    endpoints_source_checked = FALSE, # NOT established at individual trial level
    overlap_source_checked = FALSE,
    c1_c2_compared = FALSE, # original supplement has not been retrieved/inspected
    source_note = ifelse(is_barer,
      "Original review reports multi-arm Barer comparisons with a common control; exclude from independent-row illustration",
      "Package counts compatible only; participant identity and timepoints not trial-verified"),
    stringsAsFactors=FALSE
  )
}

check_supplement <- function(d, folder) {
  f <- file.path(folder, "data", "supplement_C1_C2_transcribed.csv")
  if (!file.exists(f)) return("Original C1/C2 supplement not obtained: no numerical reconciliation performed.")
  s <- read.csv(f, stringsAsFactors=FALSE)
  keys <- c("source_row", "studyID", "st_D", "sc_D", "nt_D", "nc_D", "st_DD", "sc_DD", "nt_DD", "nc_DD")
  if (!all(keys %in% names(s)) || anyDuplicated(s$source_row))
    stop("Supplement transcript needs unique source_row and all required columns; read README.")
  j <- merge(d[,keys], s[,keys], by="source_row", suffixes=c("_pkg","_supp"))
  diffs <- character()
  for (key in keys[-c(1,2)]) {
    diffs <- c(diffs, paste(key, sum(j[[paste0(key,"_pkg")]] != j[[paste0(key,"_supp")]])))
  }
  paste("Manual supplement matches:", nrow(j), "/", nrow(d), "; differences", paste(diffs, collapse="; "),
        "; follow-up and participant-overlap still require documentation.")
}

make_real <- function(d, aud) {
  selected <- which(aud$numeric_candidate)
  stopifnot(length(selected)==9, sum(aud$included_in_primary_conditional_analysis)==7)
  m <- list()
  for (which_set in c("non_Barer_7", "all_numeric_9")) {
    ss <- if (which_set=="non_Barer_7") selected[!aud$barer_shared_control[selected]] else selected
    objects <- lapply(ss, function(i) {
      r <- d[i,]
      nested_effects(r$nt_D, r$st_D, r$st_DD, r$nc_D, r$sc_D, r$sc_DD)
    })
    err <- max(vapply(objects, function(o) max(o$checks), numeric(1)))
    if (err>1e-9) stop("Nested/transform identity failed in real-data inputs")
    yy <- do.call(rbind,lapply(objects, function(o) o$y_RR))
    for (method in METHODS) {
      SS <- lapply(objects, function(o) o$transformed$S[[method]])
      fit <- gls_fit(yy,SS)$summary
      m[[length(m)+1]] <- cbind(subset=which_set, k=length(ss), method=method, fit,
                               source_verified=FALSE, model="common_effect_GLS")
    }
  }
  ans <- do.call(rbind,m)
  rownames(ans) <- NULL
  ans
}

run_replay <- function(folder, reps=2000L) {
  K <- 20L
  replay <- load_replay(folder, reps)
  sc <- simulation_scenarios()
  metrics <- list()
  valid_error <- 0
  for (i in seq_len(nrow(sc))) {
    raw <- scenario_counts(sc[i,], reps, K, replay, 20261109L)
    corrected <- correct_joint_counts(raw$ct,raw$cc)
    for (reporting in c("all_OR","mixed_OR_RR")) {
      q <- batch_quantities(corrected$ct,corrected$cc,reporting)
      valid_error <- max(valid_error, validate_batch(q, corrected$ct, corrected$cc))
      f <- lapply(METHODS, function(m) batch_gls(q$y1,q$y2,q$S[[m]]))
      names(f) <- METHODS
      truth <- c(log(.8), log(.9), log(.8)-log(.9))
      for (m in METHODS) {
        result <- simulation_metrics(f[[m]], f$full, truth, sc[i,], reporting, m, reps)
        metrics[[length(metrics)+1L]] <- result
      }
    }
    cat(sprintf("Replayed scenario %d/4: %s, n=%d, reps=%d\n",i,sc$regime[i],sc$n_per_arm[i],reps))
  }
  if (valid_error>1e-9) stop("Simulation batch vs scalar validation failed")
  x <- do.call(rbind,metrics)
  x <- x[,c("scenario","regime","n_per_arm","reporting","method","target","reps","bias","rmse",
           "mean_reported_se","coverage","coverage_mcse","mean_interval_width","failure_rate")]
  rownames(x)<-NULL
  x
}

run_main <- function() {
  args<-commandArgs(trailingOnly=TRUE)
  skip <- "--skip-sim" %in% args
  loc<-which(args=="--reps")
  reps<-if(length(loc)) as.integer(args[loc[1]+1L]) else 2000L
  if(is.na(reps)||reps<2L||reps>2000L) stop("--reps must be 2..2000; replay data only")
  argv<-commandArgs(trailingOnly=FALSE)
  f<-grep("^--file=",argv,value=TRUE)
  base<-if(length(f)) dirname(normalizePath(sub("^--file=","",f[1]))) else getwd()
  dest<-file.path(base,"results")
  dir.create(dest,recursive=TRUE,showWarnings=FALSE)
  # prevent stale full simulation results surviving a skip run
  for (name in c("source_audit.csv","pooled_comparison.csv","simulation_summary.csv","run_status.txt")) {
    path<-file.path(dest,name)
    if(file.exists(path)) unlink(path)
  }
  d<-read.csv(file.path(base,"data","geeganage_package.csv"),stringsAsFactors=FALSE)
  expected<-"a311056fce73bad8753e0c46fe9c5e1a"
  if(nrow(d)!=17L || unname(tools::md5sum(file.path(base,"data","geeganage_package.csv")))!=expected)
    stop("Package dataset changed; inspect provenance before rerunning")
  aud<-make_audit(d)
  write.csv(aud,file.path(dest,"source_audit.csv"),row.names=FALSE)
  real<-make_real(d,aud)
  write.csv(real,file.path(dest,"pooled_comparison.csv"),row.names=FALSE)
  info <- c("R execution completed; results are conditional and not a clinical meta-analysis.",
    "Package dictionary: D = early death <=1 month; DD = early death or deterioration <=1 month.",
    "Package overview uses DD=death or disability; endpoint mapping is not fully reconciled.",
    "Cochrane review separately defines early death/deterioration and end-of-trial death/dependency.",
    "Source: https://pmc.ncbi.nlm.nih.gov/articles/PMC7120409/",
    "Wei-Higgins paper Section 6 and Table C1/C2 were not numerically reconciled.",
    "Numerical candidates: 9/17; primary illustrative subset excludes two Barer rows: 7/17.",
    "Barer multi-arm trials have shared controls; all_numeric_9 treats rows as independent (not justified).",
    "Non-Barer independence and participant overlap are also unverified.",
    check_supplement(d,base))
  if (!skip) {
    sim<-run_replay(base,reps)
    write.csv(sim,file.path(dest,"simulation_summary.csv"),row.names=FALSE)
    info<-c(info,sprintf("Simulation replay: %d distinct meta-analysis datasets across four scenarios; two reporting patterns reuse each dataset.",4*reps),
             "Zero between-study heterogeneity; Wald normal intervals and specified continuity correction.")
  } else info<-c(info,"Simulation skipped by request; no stale simulation_summary.csv left in results.")
  writeLines(info,file.path(dest,"run_status.txt"))
  cat("Completed. Files in results/: source_audit.csv, pooled_comparison.csv, run_status.txt",
      if(!skip) ", simulation_summary.csv" else "", "\n",sep="")
}

if (sys.nframe()==0) run_main()
