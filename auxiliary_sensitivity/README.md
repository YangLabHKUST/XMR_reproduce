# Figure 5: auxiliary-dataset sensitivity

Hold the BBJ target GWAS fixed and compare FinnGen and Pan-UKB auxiliary exposure datasets.
Supply the 68-trait manifest, formatted GWAS and shared LD resources in [DATA.md](../DATA.md).

## Analysis

```bash
Rscript auxiliary_sensitivity/code/run_analysis.R
```

The entry point runs these three method scripts, which may also be run individually:

```bash
Rscript auxiliary_sensitivity/code/run_XMR_Fin_PanUKB_to_BBJ.R
Rscript auxiliary_sensitivity/code/run_CAMeRa_Fin_PanUKB_to_BBJ.R
Rscript auxiliary_sensitivity/code/run_MRBEE_TL_Fin_PanUKB_to_BBJ.R
```

XMR estimates its background for each trait pair with zero cross-population intercepts,
then selects IVs and fits the model. Backgrounds and summaries are under
`results/crossPopIntercept0/`.

The MRBEE/TL runner estimates source and target background intercept matrices when the
corresponding cache is absent. It saves them to `results/background_MRBEE_TL_BBJ/` and
produces both MRBEE and MRBEE-TL estimates. CAMERA and MRBEE/TL summaries are under
`results/masterClump5e-05/`. No precomputed background file is required to start the analysis.

## Plotting

```bash
Rscript auxiliary_sensitivity/plot/plot_auxiliary_sensitivity_four_methods_two_thresholds.R
```

The analysis evaluates four IV thresholds. The figure uses 5e-5 and 5e-8 and compares the
same 168 exposure–outcome pairs (three exposures and 56 outcomes). The reference set is
selected using XMR availability across both auxiliary datasets and all four thresholds;
retain the full four-threshold XMR table when preparing plot inputs. Figures go to `figures/`.
