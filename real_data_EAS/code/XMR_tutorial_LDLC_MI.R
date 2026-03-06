# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .R
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.15.2
#   kernelspec:
#     display_name: R
#     language: R
#     name: ir
# ---

# # A tutorial of MR analysis with XMR 
#
# We illustrate the XMR pipeline through a case study examining the causal effect of LDL cholesterol (LDLC) on myocardial infarction (MI) in East Asians (EAS), with Europeans (EUR) serving as the auxiliary population.

# ## step 0: Install and load XMR package

# install.packages("devtools")
devtools::install_github("YangLabHKUST/XMR")

library(XMR)

# ## step1: Preprocess and format GWAS data

# ### Data Sources
#
# Referring to the [BBJ website](https://pheweb.jp/downloads), we download the corresponding GWAS data of LDL-C and MI for East Asians and Europeans:
#
# **LDL-C:**
# - EAS (BBJ): [hum0197.v3.BBJ.LDLC.v1.zip](https://humandbs.dbcls.jp/files/hum0197/hum0197.v3.BBJ.LDLC.v1.zip)
# - EUR (UKB): [Neale Lab UK Biobank](http://www.nealelab.is/uk-biobank); [UKB Results manifest](https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/edit?gid=178908679#gid=178908679), phecode `30780`
#
# **MI:**
# - EAS (BBJ): [hum0197.v3.BBJ.MI.v1.zip](https://humandbs.dbcls.jp/files/hum0197/hum0197.v3.BBJ.MI.v1.zip)
# - EUR (UKB+FG): [hum0197.v3.EUR.MI.v1.zip](https://humandbs.dbcls.jp/files/hum0197/hum0197.v3.EUR.MI.v1.zip)

# Our analysis requires properly formatted GWAS summary statistics for both the exposure (LDL-C) and the outcome (MI), in both EAS and EUR populations.
#
# The `format_data()` function from the **MRAPSS** package standardises raw GWAS files into a uniform format. It requires, at minimum:
#
# - **SNP identifier** (rs number)
# - **Effect and non-effect alleles**
# - **A signed summary statistic** (e.g., effect size β, or odds ratio) together with its standard error or *p*-value, from which a *z*-score can be derived
# - **Sample size**
#
# Optionally, supplying **allele frequency** (`freq_col`) and **imputation quality** (`info_col`) enables additional QC filtering. 
#
# The function automatically performs several quality-control steps, including:
#
# 1. Retaining only HapMap 3 SNPs
# 2. Removing SNPs with minor allele frequency below a user-specified threshold
# 3. Excluding non-standard or ambiguous alleles (e.g., G/C, A/T)
# 4. Excluding the MHC region (chr6: 26–34 Mb)
# 5. Removing SNPs whose χ² exceeds a data-adaptive cap (default: max(*N*/1000, 80))
#
# After formatting, each dataset will contain the columns: `SNP`, `A1`, `A2`, `Z`, `chi2`, `P`, and `N`.

# Key parameter choices in `format_data()`:
#
# | Parameter | Meaning |
# |-----------|---------|
# | `snp_col` | Column containing rs IDs |
# | `b_col` / `se_col` | Effect size and its standard error |
# | `A1_col` / `A2_col` | Effect allele / non-effect allele |
# | `p_col` | *p*-value |
# | `freq_col` | Allele frequency |
# | `n` / `n_col` | Total sample size / column containing sample size |
# | `min_freq` | Remove SNPs with MAF below the threshold |

# #### Assigning rs IDs with `variants.tsv`
#
# Some GWAS summary statistics use genomic-position-based variant identifiers (e.g., `chr:pos:ref:alt`) rather than rs numbers. 
#
# Since `format_data()` filters to HapMap 3 SNPs by rs ID, we first need to map these identifiers to rs IDs:
#
# `variants.tsv` — used for **LDL-C (EUR)** from the Neale Lab UK Biobank GWAS.
#
# This annotation file accompanies the Neale Lab Round 2 results; its rows are aligned one-to-one with the GWAS file, making the mapping straightforward. 
#
# It can be downloaded from the link at the top of the [UKB Results manifest](https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/edit?gid=178908679#gid=178908679):
#
# ```bash
# wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/annotations/variants.tsv.bgz \
#      -O variants.tsv.bgz
# gunzip variants.tsv.bgz
# ```

setwd('./real_data_EAS') # follow the folder structure for reproduction to access data

# +
library(data.table)
library(MRAPSS)

### LDLC (BBJ)
exp_BBJ_address = "./raw_data/BBJ/hum0197.v3.BBJ.LDLC.v1/GWASsummary_LDLC_Japanese_SakaueKanai2020.auto.txt.gz" # replace with your own path
exp_BBJ = fread(exp_BBJ_address)
exp_BBJ = format_data(exp_BBJ,
                      snp_col = "SNP",
                      b_col = "BETA",
                      se_col = "SE",
                      A1_col = "ALLELE1", # effect allele
                      A2_col = "ALLELE0", # non-effect allele
                      p_col = "P_BOLT_LMM_INF", # 
                      freq_col = "A1FREQ",
                      n = 72866, # sample size
                      min_freq = 0.01)


### MI (BBJ)
out_BBJ_address = "./raw_data/BBJ/hum0197.v3.BBJ.MI.v1/GWASsummary_MI_Japanese_SakaueKanai2020.auto.txt.gz"  # replace with your own path
out_BBJ = data.table::fread(out_BBJ_address)
out_BBJ = format_data(out_BBJ,
                      snp_col = "SNPID",
                      b_col = "BETA",
                      se_col = "SE",
                      A1_col = "Allele2", # effect allele
                      A2_col = "Allele1", # non-effect allele
                      p_col = "p.value",
                      freq_col = "AF_Allele2",
                      n = 161206,
                      min_freq = 0.01)


### LDLC (EUR)
variants = fread('./raw_data/EUR/variants.tsv')
exp_EUR_address = "./raw_data/EUR/30780.tsv"  # replace with your own path
exp_EUR = fread(exp_EUR_address)

check1 = sum(exp_EUR$variant != variants$variant)
if(check1 != 0){
    print('Need more careful alignment with variants.tsv... skip')
    next
}

exp_EUR$SNP = variants$rsid
exp_EUR$A1 = variants$alt # effect allele
exp_EUR$A2 = variants$ref # non-effect allele
exp_EUR = format_data(exp_EUR,
                      snp_col = "SNP",
                      b_col = "beta",
                      se_col = "se",
                      A1_col = "A1",
                      A2_col = "A2",
                      p_col = "pval",
                      freq_col = "minor_AF",
                      n_col = 'n_complete_samples',
                      min_freq = 0.01)
# -

write.table(exp_BBJ, file="./BBJ_30780", quote = F, col.names = T, row.names = F)
write.table(out_BBJ, file="./BBJ_411.2", quote = F, col.names = T, row.names = F)
write.table(exp_EUR, file="./EUR_30780", quote = F, col.names = T, row.names = F)

# ## step2: Estimate background parameters
#
# Both Omega and C can be estimated via bivariate LD score regression. Here we use the functions from the `XMAP` package. The population-specific LD score file `LDscore_eas_brit_afr.txt` is available in the data provided in the `Reproducibility` section.

# #### `estimate_bgparas()`: Estimate background parameters via two-stage bivariate LDSC
#
# The function `estimate_bgparas()` estimates background parameters using a two-stage bivariate LD score regression. It is included in the `XMR` package and can be called directly after loading the library.
#
# #### Parameters
#
# | Parameter | Description | Default |
# |---|---|---|
# | `ldscore0` | A `data.frame` of LD scores (see format requirements below). | — |
# | `exposure_code` | Trait code for the exposure, e.g. `"30780"`. | — |
# | `outcome_code` | Trait code for the outcome, e.g. `"411.2"`. | — |
# | `pop_ref` | Label for the auxiliary (reference) population, e.g. `"EUR"`. | — |
# | `pop_target` | Label for the target population, e.g. `"BBJ"`. | — |
# | `gwas_dir` | Directory containing the formatted GWAS summary statistics files. | — |
# | `snp_col` | Column name for SNP rsID in `ldscore0`. | `"rsid"` |
# | `a1_col` | Column name for effect allele in `ldscore0`. | `"allele1"` |
# | `z2_threshold` | Z² threshold for filtering SNPs with large effects in stage 1. | `30` |
#
# #### Input format
#
# **GWAS files** should follow this convention:
#
# - File names take the form `{pop}_{trait_code}`, e.g. `EUR_30780`, `BBJ_411.2`.
# - Each file must contain at least the following columns: `SNP`, `Z`, `N`, `A1`.
#
# **LD score `data.frame`** (`ldscore0`) must contain the following columns, where the column names must be consistent with the `pop_ref` and `pop_target` labels you pass to the function:
#
# - A SNP identifier column (default name `rsid`).
# - An effect allele column (default name `allele1`).
# - Within-population LD score columns named after each population label, e.g. `EUR` and `BBJ`.
# - Cross-population LD score columns named `{popA}_{popB}`, including all four combinations: `EUR_EUR`, `EUR_BBJ`, `BBJ_EUR`, and `BBJ_BBJ`.
#
# #### Return value
#
# A `data.frame` with 6 rows — one for each unique (trait, population) pair — containing the following columns: `trait1`, `trait2`, `pop1`, `pop2`, `omega11`, `omega12`, `omega22`, `c1`, `c12`, `c2`.

# +
library(XMAP)
library(data.table)
library(Matrix)

set.seed(1)

# read LD scores
ldscore <- fread("./LDscore_eas_brit_afr.txt")

# pre-process: remove ambiguous SNPs
idx_amb <- which(ldscore$allele1 == comple(ldscore$allele2))
ldscore <- ldscore[-idx_amb,]

ldscore$BBJ = ldscore$EAS
ldscore$BBJ_EUR = ldscore$EUR_EAS
ldscore$EUR_BBJ = ldscore$EUR_EAS
ldscore$EUR_EUR = ldscore$EUR
ldscore$BBJ_BBJ = ldscore$EAS

ldscore = as.data.frame(ldscore)
ldscore0 = ldscore
# -

head(ldscore0)

# +
library(XMAP)
library(data.table)

set.seed(1)


bgparas <- estimate_bgparas(
    ldscore0      = ldscore,
    exposure_code = "30780", 
    outcome_code  = "411.2",      
    pop_ref       = "EUR",         
    pop_target    = "BBJ",      
    gwas_dir      = "." 
)

print(bgparas)
# -

# ## step3: IV selection; LD clumping; Run XMR

library(reshape2)
library(dplyr)
library(magrittr)
library(Matrix)

# #### `qc_bgparas()`: Quality control of background parameter matrices
#
# After estimating the background parameters with `estimate_bgparas()`, we use
# `qc_bgparas()` to reshape the output into symmetric Omega and C matrices and
# perform several quality-control checks before passing them to `fit_XMR()`.
#
# The function carries out three checks:
#
# 1. **Negative variance** — If any diagonal entry of Omega is negative, a
#    warning is issued, indicating that the heritability estimate for the
#    corresponding trait–population combination is problematic.
# 2. **Positive semi-definiteness of C** — The intercept matrix C is tested via
#    Cholesky decomposition. A warning is issued if C is not positive
#    semi-definite.
# 3. **Correlation truncation** — Off-diagonal entries of Omega that imply a
#    correlation greater than 1 in absolute value are truncated to
#    `truncation_corr` (default 0.95) times the maximum possible value, ensuring
#    the matrix remains well-behaved.
#
# #### Parameters
#
# | Parameter | Description | Default |
# |---|---|---|
# | `bgparas` | A `data.frame` returned by `estimate_bgparas()`. | — |
# | `exposure_code` | Trait code for the exposure, e.g. `"30780"`. | — |
# | `outcome_code` | Trait code for the outcome, e.g. `"411.2"`. | — |
# | `pop_ref` | Reference population label, e.g. `"EUR"`. | — |
# | `pop_target` | Target population label, e.g. `"BBJ"`. | — |
# | `truncation_corr` | Cap for off-diagonal correlations exceeding 1 in absolute value. | `0.95` |
#
# #### Return value
#
# A list with two elements:
#
# - **`Omega0`** — A 3 × 3 named symmetric matrix of background genetic covariance estimates.
# - **`C`** — A 3 × 3 named symmetric matrix of LDSC intercept estimates.

# +
colnames(bgparas) = c('trait1','trait2','pop1','pop2','Omega11','Omega12','Omega22','c1','c12','c2')

bg <- qc_bgparas(
  bgparas       = bgparas,
  exposure_code = "30780",
  outcome_code  = "411.2",
  pop_ref       = "EUR",
  pop_target    = "BBJ"
)

Omega0 <- bg$Omega0
C      <- bg$C
# -

# ### Allele harmonization with `align_data()`
#
# GWAS summary statistics from different sources may code effect alleles differently. Before analysis, we must ensure that effect sizes across all datasets refer to the same allele. 
#
# The `align_data()` function harmonizes two datasets at a time by flipping effect sizes where necessary so that allele directions are consistent.
#
# We use the target-population exposure GWAS as the anchor and call `align_data()` twice — once to harmonize the outcome GWAS and once to harmonize the auxiliary-population exposure GWAS — so that all three datasets end up on a common allele orientation.

# +
exposure_code <- "30780"
outcome_code  <- "411.2"

# p-value threshold for IV selection (applied to the auxiliary population)
threshold <- 5e-5

#  Rename columns to match the format expected by XMR functions:
#    L2.pop1 = within-population LD score of auxiliary pop (EUR)
#    L2.pop2 = within-population LD score of target pop (BBJ)
#    L12     = cross-population LD score
ldscore <- ldscore0[, c("rsid", "BBJ_BBJ", "EUR_EUR", "BBJ_EUR")]
colnames(ldscore) <- c("SNP", "L2.pop2", "L2.pop1", "L12")
ldscore0 <- ldscore

#  Load GWAS summary statistics
#  Each file should contain at least: SNP, Z, N, A1
exp.pop1 <- fread(paste0("./EUR_", exposure_code))   # exposure GWAS in auxiliary pop (EUR)
exp.pop2 <- fread(paste0("./BBJ_", exposure_code))   # exposure GWAS in target pop (BBJ)
out.pop2 <- fread(paste0("./BBJ_", outcome_code))    # outcome GWAS in target pop (BBJ)

# median sample sizes for the exposure in each population
N1 <- median(exp.pop1$N)   # auxiliary pop sample size
N2 <- median(exp.pop2$N)   # target pop sample size

# Harmonize effect allele directions across GWAS datasets

paras1 <- align_data(dat1 = exp.pop2,
                     dat2 = out.pop2,
                     trait1.name = paste0("BBJ_", exposure_code),
                     trait2.name = paste0("BBJ_", outcome_code),
                     h2.fix.intercept = FALSE,
                     LDSC = FALSE)

paras2 <- align_data(dat1 = exp.pop2,
                     dat2 = exp.pop1,
                     trait1.name = paste0("BBJ_", exposure_code),
                     trait2.name = paste0("EUR_", exposure_code),
                     h2.fix.intercept = FALSE,
                     LDSC = FALSE)

# Rename harmonized columns for clarity
dat1 <- paras1$dat
colnames(dat1) <- c("SNP", "A1", "A2",
                     "b.exp.pop2", "b.out.pop2",
                     "se.exp.pop2", "se.out.pop2",
                     "pval.exp.pop2", "pval.out.pop2")

dat2 <- paras2$dat
colnames(dat2) <- c("SNP", "A1", "A2",
                     "b.exp.pop2", "b.exp.pop1",
                     "se.exp.pop2", "se.exp.pop1",
                     "pval.exp.pop2", "pval.exp.pop1")

# Merge all three GWAS datasets and LD scores into a single data.frame
dat <- merge(dat1, dat2[, c("SNP", "A1", "A2",
                             "b.exp.pop1", "se.exp.pop1", "pval.exp.pop1")])
dat <- merge(dat, ldscore0, by = "SNP")
# -

# ### LD clumping in the target population
#
# To obtain approximately independent instrumental variables, we perform LD clumping using PLINK with the target-population LD reference panel. 
#
# SNPs are first filtered by the IV selection threshold (applied to the auxiliary-population p-values), and then clumped at an r² threshold of 0.001 within a 1 Mb window using target-population LD structure.
#
# Both the PLINK and the 1000 Genomes LD reference panels are provided in the `Reproducibility` section.
# If you need to install PLINK yourself, pre-compiled binaries for Linux, macOS, and Windows are available at [https://www.cog-genomics.org/plink/1.9/](https://www.cog-genomics.org/plink/1.9/). Simply download the appropriate archive, extract the binary, and supply its path to the `plink_bin` argument.

# +
#  LD clumping is performed in the target population (EAS/BBJ) to prune SNPs that are in linkage disequilibrium, retaining approximately independent instruments.
#
#  - IV pre-selection: keep SNPs with p < threshold in the auxiliary pop
#  - Clumping window: 1000 kb
#  - Clumping r²:     0.001
#  - bfile:           target-population LD reference panel (EAS)
#  - plink_bin:       path to the PLINK


clumped.MR <- MRAPSS::clump(subset(dat, pval.exp.pop1 < threshold),
                           IV.Threshold = 1,
                           SNP_col      = "SNP",
                           pval_col     = "pval.exp.pop2",
                           clump_kb     = 1000,
                           clump_r2     = 0.001,
                           bfile        = "../nc/1kg_pops/EAS_phase3_hg38_autosomes",
                           plink_bin    = "../nc/plink_Linux"
                          )
clumped.MR$Threshold <- threshold
# -

# ### Modified threshold for selection-bias correction
#
# LD clumping introduces additional selection bias.
# And in XMR, IV selection is based on a p-value threshold applied to the auxiliary population (e.g. EUR), while LD clumping is performed in the target population (e.g. BBJ).
#
# The modified threshold accounts for the genetic relationship between the two populations — reflected in their relative GWAS sample sizes — and adjusts the nominal threshold accordingly:
# This corrected threshold is then passed to `fit_XMR()` so that the likelihood properly reflects the cross-population selection process.

# +
# Compute modified threshold
#  Because LD clumping is applied to the target population while IV selection is based on a p-value threshold in the auxiliary population, 
#  the modified threshold used to address additional selection bias from LD clumping accounts for the relationship between the two populations (reflected in their relative sample sizes N1 and N2)

t0 <- abs(qnorm(threshold / 2))
dt <- 0.13 / sqrt(N2 / N1)
modified_threshold <- 2 * (1 - pnorm(abs(t0 + dt)))
# -

# ### Fitting the XMR model
#
# With the clumped instruments, background parameter matrices (`Omega0` and `C`), and the modified threshold in hand, we fit the XMR model.
#
# When the number of IVs is small (fewer than 100), tighter convergence tolerances are used for more precise estimation; otherwise the tolerances are relaxed slightly.
#
# The returned object contains the causal effect estimate (`beta`), its standard error (`beta.se`), the corresponding p-value (`beta.pvalue`), and additional diagnostics such as the number of valid instruments (`nvalid`).

# +
# Fit XMR model
#  When the number of IVs is small (< 100) we use tighter convergence tolerances (1e-08) for more precise estimation; 
#  otherwise we relax them slightly (1e-07) to speed up computation.

if (nrow(clumped.MR) < 100) {
  XMR_res <- fit_XMR(data      = clumped.MR,
                     C         = C,
                     Omega0    = Omega0,
                     Threshold = modified_threshold,
                     tol1      = 1e-08,
                     tol2      = 1e-08,
                     min_thres = 1e-2)
} else {
  XMR_res <- fit_XMR(data      = clumped.MR,
                     C         = C,
                     Omega0    = Omega0,
                     Threshold = modified_threshold,
                     tol1      = 1e-07,
                     tol2      = 1e-07,
                     min_thres = 1e-2)
}

#  Collect results

res_xmr <- data.frame(exposure            = exposure_code,
                      outcome             = outcome_code,
                      beta                = XMR_res$beta,
                      beta.se             = XMR_res$beta.se,
                      beta.pvalue         = XMR_res$beta.pvalue,
                      tau.sq              = XMR_res$tau.sq,
                      sigmax11            = XMR_res$SigmaX[1, 1],
                      sigmax12            = XMR_res$SigmaX[1, 2],
                      sigmax22            = XMR_res$SigmaX[2, 2],
                      nIV                 = XMR_res$nIV,
                      IVsignal.sum        = XMR_res$post$IVsignal.sum,
                      nvalid              = XMR_res$nvalid,
                      pi0                 = XMR_res$pi0,
                      Threshold           = threshold,
                      exp_N1              = N1,
                      exp_N2              = N2,
                      modified_threshold  = modified_threshold
                    )
# -

res_xmr

# Our results demonstrate that in the EAS population (BBJ), LDL cholesterol (LDLC) exhibits a significant positive causal effect on myocardial infarction (MI), with a 
# highly significant p-value of $3.799 \times 10^{-7}$, highlighting the strong statistical power of our proposed method. 
# The discovery also aligns with established literature. 
#
# Furthermore, the model fitting procedure can be completed in around one minute on a single Intel® Xeon® Gold 6152 CPU @ 2.10 GHz, demonstrating the computational efficiency of our approach.
