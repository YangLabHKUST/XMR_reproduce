# Figure 7: BBJ–TPMI

Use EUR auxiliary exposure data with BBJ and TPMI targets. The formatted TPMI files use the
`TW_` prefix. The Pan-UKB source-outcome files used by comparison methods are specified in
[DATA.md](../DATA.md).

## Analysis

```bash
Rscript BBJ_TPMI/code/run_analysis.R
```

The entry point runs these method scripts in order. They can also be run individually:

```bash
Rscript BBJ_TPMI/code/run_XMR_EAS_BBJ_TPMI.R
Rscript BBJ_TPMI/code/run_CAMERA_EUR_to_BBJ_TPMI.R
Rscript BBJ_TPMI/code/run_MRBEE_TL_EUR_to_BBJ_TPMI.R
Rscript BBJ_TPMI/code/run_MRAPSS_BBJ_TPMI.R
Rscript BBJ_TPMI/code/run_MRAPSS_EUR_PanUKB_5e-05.R
```

XMR estimates a genome-wide three-trait background with zero cross-population intercepts.
The MR-APSS runners estimate their population-specific backgrounds. MRBEE/TL estimates
source and target intercept matrices when its cache under
`results/background_MRBEE_TL_EUR_PanUKB_to_BBJ_TPMI/` is absent.
No precomputed background file is required to start these analyses.

`shared_target_IVs_BBJ_TPMI.R` is a dependency of the target MR-APSS runner.
XMR, CAMERA and MRBEE/TL write to `results/masterClump5e-05/`; target MR-APSS writes to
`results/sharedIV_masterClump5e-05/MRAPSS/`; EUR MR-APSS writes to `results/EUR/`.
When using previously generated tables, place them in this layout. Keep both thresholds in
the XMR `two_thresholds` tables; plotting selects 5e-5 where required.

## Plotting

```bash
Rscript BBJ_TPMI/plot/plot_BBJ_TPMI_UKBB_MR_comparison_ABCDE.R
```

This entry point automatically combines panels A–E using the other R scripts in `plot/`.
It reads result tables directly and saves PDF, PNG and plot-data tables in `figures/`.

Panel E uses each method's own available fits. In TPMI, XMR, CAMERA, MRBEE and MRBEE-TL
use 105 pairs, and MR-APSS uses its 96 available pairs. Failed MR-APSS fits do not remove
pairs from the other methods. Bonferroni correction uses the 105 planned pairs.
