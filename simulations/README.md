# Figure 3: simulations

`code/baseline/` and `code/reduced_confounding/` contain separate numerical workflows:
LD/SNP preparation, effect generation, summary-statistic generation, background estimation,
and analysis with XMR, CAMERA, MRBEE/MRBEE-TL, TEMR and MR-APSS.

## Prepare LD inputs

Place the White British and EAS ViLMA reference blocks at the locations in [DATA.md](../DATA.md).
From the repository root, create the shared chromosome 2 and 19 SNP and LD-score tables:

```bash
python3 simulations/code/baseline/step1_extract_vilma_snps_bychr.py --ld-dir simulations/data/vilma_ld --chromosomes 2 19
python3 simulations/code/baseline/step2_compute_vilma_ldscores.py --ld-dir simulations/data/vilma_ld --chromosomes 2 19
```

These prepared inputs are shared by both simulation profiles. When using a custom data
location, pass that LD directory to both preparation commands and set `XMR_SIM_DATA_ROOT`
to its parent directory for the analysis launcher.

## Analysis

[fig3_settings.tsv](fig3_settings.tsv) lists the 15 configurations needed for the main figure:

| Profile | Target sample size | Auxiliary sample sizes | Exposure-effect correlations |
| --- | --- | --- | --- |
| Baseline | 6,000 | 20,000; 40,000; 60,000 | 0.2; 0.5; 0.8 |
| Baseline | 8,000 | 40,000 | 0.2; 0.5; 0.8 |
| Reduced confounding | 6,000 | 40,000 | 0.2; 0.5; 0.8 |

Each configuration evaluates causal effects 0, 0.1, 0.2, 0.3 and 0.4 with 50 replicates,
IV thresholds 5e-5 and 5e-8, and master seed 1234. The separate profile parameter files
specify the background and sampling-error correlations.

Inspect the complete command list without running simulations:

```bash
bash simulations/code/run_fig3.sh --all --dry-run
```

Run all configurations, or one configuration:

```bash
bash simulations/code/run_fig3.sh --all
# Arguments: profile, target sample size, auxiliary sample size, correlation.
bash simulations/code/run_fig3.sh baseline 6000 40000 0.8
```

The launcher runs data generation, background estimation and method fitting in that order.
All stages use the same configuration output root. Method runners receive the background
path explicitly; precomputed background parameters are not required.

Baseline outputs are under
`results/baseline/naux<auxiliary_n>_ntarget<target_n>_rho<rho_label>/`, where `rho_label`
is `0p2`, `0p5` or `0p8`. Reduced-confounding outputs use the corresponding profile directory
and append `_bgRhoLinked_uxvy0.1_vyux0.1_easCorr0.2_C1.05_0.05_1.05` to the configuration name.
Each configuration contains `true_effects/`, `summary_statistics/`, background directories
and `method_results/`.

## Result summaries and plotting

For a general summary of one configuration:

```bash
Rscript simulations/plot/baseline/evaluate_simulation_study.R --root simulations/results/baseline/naux40000_ntarget6000_rho0p8
```

Use the reduced-confounding evaluation script for that profile and pass its complete output
root, including the suffix above. Evaluation writes `evaluation/` and `figures/` inside the
configuration output root.

The final Figure 3 composition script is pending. The included evaluation scripts do not
implement its pooled MR-APSS reference curves or replicate-cluster bootstrap intervals.
