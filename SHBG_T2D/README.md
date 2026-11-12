# Figure 6: SHBG–T2D

EUR is the auxiliary population; EAS, AFR and CSA are the targets. Supply the formatted
SHBG (`30830`) and T2D (`250.2`) GWAS and the shared LD resources in [DATA.md](../DATA.md).

## Analysis

```bash
Rscript SHBG_T2D/code/run_analysis.R
```

This runs the following steps in order. They can also be run individually:

```bash
Rscript SHBG_T2D/code/prepare_background.R
Rscript SHBG_T2D/code/run_SHBG_T2D_XMR.R --arms stage2_C0
Rscript SHBG_T2D/code/run_SHBG_T2D_CAMeRa.R
Rscript SHBG_T2D/code/run_SHBG_T2D_MRAPSS.R
Rscript SHBG_T2D/code/run_SHBG_T2D_MRBEE_TL.R
Rscript SHBG_T2D/code/run_SHBG_T2D_TEMR.R
```

| Background | Computation and use | Output relative to this folder |
| --- | --- | --- |
| MR-APSS and MRBEE/TL | `prepare_background.R` estimates backgrounds for EUR, EAS, AFR and CSA. It explicitly uses `cross_intercept_zero = FALSE` and applies the covariance QC used by these analyses. The two method runners read the resulting population-indexed `C` and `Omega0` lists. | `results/SHBG-T2D-bgparas.RData` |
| XMR | The XMR runner estimates backgrounds with `cross_intercept_zero = TRUE`, selects IVs and fits the `stage2_C0` analysis. | `results/xmr_pipeline/` |

The XMR summary for plotting is written to
`results/crossPopIntercept0/SHBG_T2D_XMR_EUR_selected_FEMA_P_target_clump_noEffectHetQC.tsv`.
The comparison methods write to their respective directories under `results/`.

## Plotting

```bash
Rscript SHBG_T2D/plot/plot_SHBG_T2D_selected_methods_forest_by_method.R
```

The plot reads the method summary tables and CAMERA's
`results/CAMeRa_results_clump10000kb_V2/threshold_5e-08/FEMA_IV/CAMeRa_harmonised_SHBG_T2D.tsv`
for IV counts. Keep that table with the CAMERA results. Figures are written to `figures/`.
