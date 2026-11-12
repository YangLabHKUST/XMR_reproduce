# Figure 3: simulation plots

These R scripts create the twelve-panel simulation figure from completed method fits.
They summarize results and draw the figure; data generation and method fitting are run
separately with `simulations/code/run_fig3.sh`.

## Plotting

Run the commands below from the repository root. The input directory must contain `baseline/`
and `reduced_confounding/`, with the `simulation_*_results.tsv` files produced by the
15 configurations in [fig3_settings.tsv](../fig3_settings.tsv). Tables are found recursively, so
configuration directory names can vary. Keep one result set per configuration.

The plots use `data.table`, `ggplot2`, `cowplot` and `scales`; see the
[software dependencies](../../environment/README.md).

```bash
Rscript simulations/plot/prepare_figure3.R simulations/results simulations/results/figure3
Rscript simulations/plot/plot_figure3.R simulations/results/figure3 simulations/figures
```

The first command writes `qq.tsv`, `power_bars.tsv`, `power_curves.tsv` and
`fit_counts.tsv`. The second writes `Figure_3.pdf`, `Figure_3.png` and `sessionInfo.txt`.
The PDF is vector-based; the PNG is 300 dpi. Both use a 22 × 17.5 inch canvas.

`prepare_figure3.R` also accepts an optional third argument: a per-replicate TSV for
the auxiliary-sample-size comparison. This supports saved result exports when the
individual fourth-column result files are not available. It requires `method`,
`threshold`, `beta_true`, `rho`, `n_aux`, `n_target`, `replicate`, `status`, `beta`, `se`
and `p`. Use `NA` for the auxiliary sample size of an already deduplicated MR-APSS
reference. This argument supplies replicate-level fits, not precomputed power estimates.

## Panels and statistics

| Panels | Setting |
| --- | --- |
| A/B/C | Baseline; target N = 6,000; auxiliary N = 40,000; rho = 0.8 |
| E/F/G | Baseline; target N = 8,000; auxiliary N = 40,000; rho = 0.8 |
| I/J/K | Reduced confounding; target N = 6,000; auxiliary N = 40,000; rho = 0.8 |
| D/H/L | Baseline; target N = 6,000; auxiliary N = 20,000/40,000/60,000; rho = 0.2/0.5/0.8 |

The first two columns show null QQ plots at IV-selection thresholds 5e-5 and 5e-8;
MR-APSS is omitted from the 5e-8 QQ panels. Ranks are calculated separately within
each setting and method. Pointwise null envelopes use beta order statistics.
Observed -log10(P) values above 3 are displayed at 3; uncapped values are retained
in `qq.tsv`.

C/G/K compare XMR and MR-APSS at 5e-5; K also shows CAMERA and target-population
MRBEE at 5e-8. D/H/L compare XMR auxiliary sample sizes and MR-APSS, all at 5e-5.
The full panel mapping is in [panel_map.tsv](panel_map.tsv).

A valid fit has `status == "ok"`, finite beta, finite positive SE, and finite P in
[0, 1]. Power is the fraction of valid fits with P < 0.05; beta = 0 gives the null
rejection rate. The scripts check the expected 50 replicates per displayed setting.

MR-APSS power pools the three rho settings at fixed background and target sample size,
using one auxiliary-size copy. This gives 150 records in 50 replicate clusters.
C/G/K use 95% replicate-cluster bootstrap percentile intervals (20,000 draws; seed
20261029). Bootstrap draws use a fixed setting order, including slots for
unplotted target sizes, so supplying only the fifteen Figure 3 configurations
does not change these intervals. D/H/L and the other power series retain the exact
binomial intervals used in the figure.

Verified with R 4.4.2, data.table 1.18.2.1, ggplot2 4.0.2, cowplot 1.2.0 and scales 1.4.0.
Font rendering can vary between systems.
