# Software

The analysis environment uses R 4.3.3. Required R packages include XMR, XMAP, MRAPSS,
CAMeRa, MRBEE, MRBEEX, TEMR, data.table, Matrix, mvtnorm, expm, digest, ieugwasr,
GenomicRanges, IRanges, tibble, R.utils and parallel. Plotting additionally uses ggplot2, cowplot,
ggrepel and scales. PLINK is required for real-data clumping.

Use the XMR implementation at commit `e847e78f687ad0b3dc524624977deeea55e91d37`:

```r
remotes::install_github("YangLabHKUST/XMR@e847e78f687ad0b3dc524624977deeea55e91d37")
```

`R_packages.tsv` lists analysis-environment versions and
available source commits for the comparison methods. It records the analysis dependencies;
it is not a complete environment lockfile.

Simulation Python scripts require NumPy and pandas, in addition to the standard library.
Use an environment containing both the R and Python dependencies. The portable launcher uses
`Rscript` and `python3` from `PATH`; `SIM_PYTHON` can select a different Python executable.
The R user-library path can be set with `R_LIBS_USER` for a local installation.

The independent-SNP MRBEE-TL helper used for simulations is included alongside the simulation
runners. Real-data runners use the MRBEEX package interface.
