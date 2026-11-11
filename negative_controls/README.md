# Figure 4: negative controls

EUR is the auxiliary population; AFR and CSA are the targets. Inputs are the formatted GWAS,
37-trait `info.csv` (35 exposures and two outcomes), LD scores and panels listed in
[DATA.md](../DATA.md).

## Analysis

```bash
Rscript negative_controls/code/run_analysis.R
```

This runs the following steps in order. The commands can also be run individually:

```bash
Rscript negative_controls/code/prepare_background.R
Rscript negative_controls/code/run_NC_XMR_EUR_to_AFR_CSA.R
Rscript negative_controls/code/run_NC_CAMeRa_EUR_to_AFR_CSA.R
Rscript negative_controls/code/run_NC_CAMeRa_EUR_AFR_FEMA_EUR.R
Rscript negative_controls/code/run_NC_CAMeRa_EUR_CSA_FEMA_EUR.R
Rscript negative_controls/code/run_NC_MRBEE_TL_EUR_to_AFR_CSA.R
Rscript negative_controls/code/run_NC_TEMR_EUR_to_AFR_CSA.R
Rscript negative_controls/code/run_NC_MRAPSS_EUR.R
```

| Background | Computation and use | Output relative to this folder |
| --- | --- | --- |
| MRBEE/TL target | `prepare_background.R` fits the six trait-pair LDSC regressions, retaining stage-1 intercepts in stage 2. The MRBEE/TL runner reads the target exposure–outcome intercept block. | `results/background_MRBEE_TL_target/<target>_<outcome>/bgparas_<exposure>` |
| MRBEE/TL EUR source | The MRBEE/TL runner estimates the EUR exposure–outcome intercept matrix when its cache is absent. | `results/background_MRBEE_TL_EUR/` |
| XMR | The XMR runner estimates its own three-trait background with zero cross-population intercepts, then selects IVs and fits XMR. | `results/crossPopIntercept0/background_XMR/` |
| EUR MR-APSS | The MR-APSS runner estimates its EUR exposure–outcome background before fitting. | `results/MRAPSS_EUR_background/` |

The preparation step writes one complete, headerless 10-column table per target/exposure/outcome
combination (140 tables). It replaces a table only after all six regressions succeed.
These comparator tables are distinct from the XMR background objects.

## Plotting

```bash
Rscript negative_controls/plot/plot_negative_control_QQ_by_method_CSA_AFR_EUR_3x5.R
```

The plot reads seven summary tables: the XMR table under `results/crossPopIntercept0/` and
six comparator tables under `results/`. It writes the Figure 4 QQ layout to `figures/`.
