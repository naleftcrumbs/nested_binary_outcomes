# Nested binary-outcome

## How to run

From the extracted folder:

```bash
Rscript run_analysis.R --reps 2000
```

For an initial check without the simulation:

```bash
Rscript run_analysis.R --skip-sim
```

## Sources

1. Wei and Higgins, *Statistics in Medicine* (2013), doi:10.1002/sim.5679.
2. Lu, *Frontiers in Psychology* (2023), doi:10.3389/fpsyg.2023.1185012.
3. Trikalinos, Hoaglin and Schmid, AHRQ Methods Research Report (2013).
4. Geeganage and Bath, *Vasoactive drugs for acute stroke*, Cochrane (2010), https://pmc.ncbi.nlm.nih.gov/articles/PMC7120409/.
5. `metavcov` package documentation: https://cran.rstudio.com/web/packages/metavcov/refman/metavcov.html.
