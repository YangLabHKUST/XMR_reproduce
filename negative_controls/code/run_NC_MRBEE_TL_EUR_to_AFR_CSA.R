#!/usr/bin/env Rscript

# Resolve repository-relative inputs independently of the launch directory.
.xmr_args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(.xmr_args)) stop("Run this entry point with Rscript.")
.xmr_root <- dirname(normalizePath(sub("^--file=", "", .xmr_args[[1L]])))
while (!file.exists(file.path(.xmr_root, ".xmr-reproduce-root"))) {
  .xmr_parent <- dirname(.xmr_root)
  if (identical(.xmr_parent, .xmr_root)) stop("Repository root not found.")
  .xmr_root <- .xmr_parent
}
setwd(.xmr_root)

############################################################
# Negative-control MRBEE-TL: EUR auxiliary exposure -> AFR/CSA targets
#
# Auxiliary population: EUR
# Target populations: AFR and CSA
# Exposures: all phenocodes in info.csv except 1747_5 and 1727
# Outcomes: 1747_5 and 1727
# IV thresholds: 5e-5 and 5e-8
# LD clumping: 10000 kb, r2 = 0.001
############################################################

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ieugwasr)
  library(MRBEE)
  library(MRBEEX)
  library(XMAP)
})

set.seed(1)

## =========================================================
## 1. Analysis parameters
## =========================================================

nc_dir <- file.path(.xmr_root, "negative_controls")
formatted_dir <- file.path(nc_dir, "formatted_data")
info_file <- file.path(nc_dir, "info.csv")
bgparas_dir <- file.path(nc_dir, "results", "background_MRBEE_TL_target")
result_dir <- file.path(.xmr_root, "negative_controls/results")
source_background_dir <- file.path(result_dir, "background_MRBEE_TL_EUR")
ldscore_file <- file.path(nc_dir, "LDscore_eas_brit_afr.txt")

# Use the EUR PLINK reference specified in DATA.md.
plink_bin <- file.path(nc_dir, "plink_Linux")
eur_bfile <- file.path(nc_dir, "1kg_pops/EUR_phase3_hg38_autosomes")

combinations <- data.table(
  auxiliary = c("EUR", "EUR"),
  target = c("AFR", "CSA")
)

thresholds <- c(5e-5, 5e-8)
clump_kb <- 10000
clump_r2 <- 0.001
gxe_p_threshold <- 5e-8
sampling_time <- 1000

result_file <- file.path(result_dir, "NC_MRBEE_TL_EUR_to_AFR_CSA_clump10000kb.tsv")

# Read the negative-control trait pairs from info.csv.
if (!file.exists(info_file)) stop("Negative-control info file not found: ", info_file)
info <- fread(info_file)
if (!"phenocode" %in% names(info)) stop("info.csv must contain phenocode")
outcomes <- c("1747_5", "1727")
exposures <- setdiff(unique(as.character(info$phenocode)), outcomes)
exposures <- exposures[!is.na(exposures) & nzchar(exposures)]
if (length(exposures) == 0L) stop("No negative-control exposures found")
if (!dir.exists(formatted_dir)) stop("Formatted-data directory not found: ", formatted_dir)
if (!file.exists(plink_bin)) stop("PLINK binary not found: ", plink_bin)
if (!file.exists(ldscore_file)) stop("LD-score file not found: ", ldscore_file)
required_eur_bfiles <- paste0(eur_bfile, c(".bed", ".bim", ".fam"))
missing_eur_bfiles <- required_eur_bfiles[!file.exists(required_eur_bfiles)]
if (length(missing_eur_bfiles) > 0L) {
  stop("Missing EUR reference files: ", paste(missing_eur_bfiles, collapse = ", "))
}

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_background_dir, recursive = TRUE, showWarnings = FALSE)
if (file.exists(result_file)) file.remove(result_file)

## =========================================================
## 2. Read formatted GWAS data
## =========================================================

find_gwas_file <- function(prefix, trait) {
  stem <- file.path(formatted_dir, paste0(prefix, "_", trait))
  files <- c(stem, paste0(stem, ".txt"), paste0(stem, ".tsv"), paste0(stem, ".gz"))
  files <- files[file.exists(files)]
  if (length(files) == 0) stop("Missing GWAS file: ", stem)
  files[1]
}

read_gwas <- function(prefix, trait) {
  dat <- fread(find_gwas_file(prefix, trait))
  required <- c("SNP", "A1", "A2", "Z", "N", "P")
  if (!all(required %in% names(dat))) {
    stop(prefix, "_", trait, " lacks columns: ", paste(setdiff(required, names(dat)), collapse = ", "))
  }

  # Construct standardized beta and SE from the formatted Z and N columns.
  dat[, `:=`(
    SNP = as.character(SNP),
    A1 = toupper(as.character(A1)),
    A2 = toupper(as.character(A2)),
    Z = as.numeric(Z),
    N = as.numeric(N),
    P = as.numeric(P)
  )]
  # Inputs are QC-filtered; construct standardized model effects.
  dat[, `:=`(beta = Z / sqrt(N), se = 1 / sqrt(N))]
  dat
}

complement_allele <- function(x) chartr("ATCG", "TAGC", x)

# Determine effect-direction flips relative to the reference alleles.
# Allele matching uses the supplied alleles without an EAF-based filter.
allele_sign <- function(ref_a1, ref_a2, a1, a2) {
  valid <- !is.na(ref_a1) & !is.na(ref_a2) & !is.na(a1) & !is.na(a2)
  same <- valid & (
    (a1 == ref_a1 & a2 == ref_a2) |
    (complement_allele(a1) == ref_a1 & complement_allele(a2) == ref_a2)
  )
  reverse <- valid & (
    (a1 == ref_a2 & a2 == ref_a1) |
    (complement_allele(a1) == ref_a2 & complement_allele(a2) == ref_a1)
  )

  sign <- rep(NA_real_, length(ref_a1))
  sign[same] <- 1
  sign[!same & reverse] <- -1
  sign
}

# FEMA uses source and target exposure associations.
harmonise_exposures <- function(source, target) {
  snps <- intersect(source$SNP, target$SNP)
  source <- source[match(snps, source$SNP)]
  target <- target[match(snps, target$SNP)]
  sign <- allele_sign(source$A1, source$A2, target$A1, target$A2)
  keep <- is.finite(sign)

  data.table(
    SNP = source$SNP[keep],
    A1 = source$A1[keep],
    A2 = source$A2[keep],
    b_source = source$beta[keep],
    se_source = source$se[keep],
    z_source = source$Z[keep],
    p_source = source$P[keep],
    b_target = target$beta[keep] * sign[keep],
    se_target = target$se[keep],
    z_target = target$Z[keep] * sign[keep],
    p_target = target$P[keep]
  )
}

# After IV selection, align outcome associations to the source exposure alleles.
add_outcomes <- function(exposure_dat, source_outcome, target_outcome, snps) {
  dat <- exposure_dat[match(snps, exposure_dat$SNP)]
  source_outcome <- source_outcome[match(snps, source_outcome$SNP)]
  target_outcome <- target_outcome[match(snps, target_outcome$SNP)]

  source_sign <- allele_sign(dat$A1, dat$A2, source_outcome$A1, source_outcome$A2)
  target_sign <- allele_sign(dat$A1, dat$A2, target_outcome$A1, target_outcome$A2)
  keep <- is.finite(source_sign) & is.finite(target_sign)

  dat <- dat[keep]
  dat[, `:=`(
    b_out_source = source_outcome$beta[keep] * source_sign[keep],
    se_out_source = source_outcome$se[keep],
    b_out_target = target_outcome$beta[keep] * target_sign[keep],
    se_out_target = target_outcome$se[keep]
  )]
  dat[complete.cases(dat)]
}

## =========================================================
## 3. Background helpers
## =========================================================

# Read target background intercepts generated by prepare_background.R.
read_saved_bgparas <- function(target, exposure, outcome) {
  filename <- file.path(
    bgparas_dir, paste0(target, "_", outcome), paste0("bgparas_", exposure)
  )
  if (!file.exists(filename)) stop("Background file not found: ", filename,
    "; run negative_controls/code/prepare_background.R first.")

  paras <- fread(filename, header = FALSE)
  if (ncol(paras) != 10L) stop("Unexpected background file format: ", filename)
  setnames(paras, c(
    "trait1", "trait2", "pop1", "pop2", "Omega11", "Omega12",
    "Omega22", "c1", "c12", "c2"
  ))

  traits <- c(
    paste0("EUR_", exposure),
    paste0(target, "_", exposure),
    paste0(target, "_", outcome)
  )
  C <- matrix(NA_real_, 3L, 3L, dimnames = list(traits, traits))
  for (i in seq_len(nrow(paras))) {
    trait1 <- as.character(paras$trait1[i])
    trait2 <- as.character(paras$trait2[i])
    if (trait1 %in% traits && trait2 %in% traits) {
      C[trait1, trait2] <- as.numeric(paras$c12[i])
      C[trait2, trait1] <- as.numeric(paras$c12[i])
    }
  }
  if (any(!is.finite(C))) stop("Background file lacks required C entries: ", filename)
  list(C = C, traits = traits, filename = filename)
}

ldscore <- fread(ldscore_file)
ldscore[, `:=`(
  rsid = as.character(rsid),
  allele1 = toupper(as.character(allele1)),
  allele2 = toupper(as.character(allele2)),
  EUR = as.numeric(EUR),
  EUR_EUR = as.numeric(EUR)
)]
ldscore <- ldscore[
  allele1 != complement_allele(allele2) & !duplicated(rsid)
]

align_z_to_ldscore <- function(dat, ld) {
  sign <- allele_sign(ld$allele1, ld$allele2, dat$A1, dat$A2)
  z <- dat$Z * sign
  z[!is.finite(sign)] <- NA_real_
  z
}

# Estimate the EUR exposure/outcome LDSC intercept matrix needed by the source
# MRBEE analysis. It cannot be taken from the three-trait target XMR background.
estimate_source_Rxy <- function(exposure_dat, outcome_dat) {
  snps <- Reduce(intersect, list(ldscore$rsid, exposure_dat$SNP, outcome_dat$SNP))
  exposure_dat <- exposure_dat[match(snps, exposure_dat$SNP)]
  outcome_dat <- outcome_dat[match(snps, outcome_dat$SNP)]
  ld <- ldscore[match(snps, ldscore$rsid)]

  z1 <- align_z_to_ldscore(exposure_dat, ld)
  z2 <- align_z_to_ldscore(outcome_dat, ld)
  valid <- is.finite(z1) & is.finite(z2) &
    is.finite(exposure_dat$N) & is.finite(outcome_dat$N) &
    is.finite(ld$EUR)
  stage1 <- valid & z1^2 < 30 & z2^2 < 30
  if (sum(stage1) < 100L) stop("Too few variants for EUR source background estimation")

  w <- 1 / pmax(ld$EUR, 1)
  fit1 <- XMAP::estimate_gc(
    data.frame(Z = z1[stage1], N = exposure_dat$N[stage1]),
    data.frame(Z = z2[stage1], N = outcome_dat$N[stage1]),
    ld$EUR[stage1], ld$EUR[stage1], ld$EUR_EUR[stage1],
    reg_w1 = w[stage1], reg_w2 = w[stage1], reg_wx = w[stage1],
    constrain_intercept = FALSE
  )
  intercept <- c(fit1$tau1$coefs[1], fit1$tau2$coefs[1], fit1$theta$coefs[1])
  if (any(!is.finite(intercept))) stop("Invalid EUR source LDSC intercept")

  fit2 <- XMAP::estimate_gc(
    data.frame(Z = z1[valid], N = exposure_dat$N[valid]),
    data.frame(Z = z2[valid], N = outcome_dat$N[valid]),
    ld$EUR[valid], ld$EUR[valid], ld$EUR_EUR[valid],
    reg_w1 = w[valid], reg_w2 = w[valid], reg_wx = w[valid],
    constrain_intercept = TRUE,
    intercept[1], intercept[2], intercept[3]
  )
  Rxy <- matrix(c(
    fit2$tau1$coefs[1], fit2$theta$coefs[1],
    fit2$theta$coefs[1], fit2$tau2$coefs[1]
  ), 2L, 2L, byrow = TRUE)
  if (any(!is.finite(Rxy)) || any(diag(Rxy) <= 0)) stop("Invalid EUR source Rxy")
  rho <- Rxy[1, 2] / sqrt(Rxy[1, 1] * Rxy[2, 2])
  if (abs(rho) >= 1) {
    Rxy[1, 2] <- Rxy[2, 1] <-
      0.95 * sign(rho) * sqrt(Rxy[1, 1] * Rxy[2, 2])
  }
  Rxy
}

get_background <- function(auxiliary, target, exposure, outcome,
                           source_exposure, source_outcome,
                           target_exposure, target_outcome) {
  saved <- read_saved_bgparas(target, exposure, outcome)
  target_exp <- paste0(target, "_", exposure)
  target_out <- paste0(target, "_", outcome)

  source_file <- file.path(
    source_background_dir,
    paste0(auxiliary, "_", exposure, "_", outcome, "_Rxy.rds")
  )
  C_source <- if (file.exists(source_file)) {
    readRDS(source_file)
  } else {
    value <- estimate_source_Rxy(source_exposure, source_outcome)
    saveRDS(value, source_file)
    value
  }

  list(
    C_source = C_source,
    C_target = saved$C[c(target_exp, target_out), c(target_exp, target_out)],
    filename = paste(source_file, saved$filename, sep = ";")
  )
}

## =========================================================
## 4. Prepare one exposure-outcome pair
## =========================================================

prepare_pair <- function(auxiliary, target, exposure, outcome) {
  # Use population-specific exposure and outcome associations for source and target.
  source_exposure <- read_gwas(auxiliary, exposure)
  source_outcome <- read_gwas(auxiliary, outcome)
  target_exposure <- read_gwas(target, exposure)
  target_outcome <- read_gwas(target, outcome)

  background <- get_background(
    auxiliary, target, exposure, outcome,
    source_exposure, source_outcome, target_exposure, target_outcome
  )

  list(
    exposure = harmonise_exposures(source_exposure, target_exposure),
    source_outcome = source_outcome,
    target_outcome = target_outcome,
    background = background,
    source_exposure_N = median(source_exposure$N, na.rm = TRUE),
    source_outcome_N = median(source_outcome$N, na.rm = TRUE),
    target_exposure_N = median(target_exposure$N, na.rm = TRUE),
    target_outcome_N = median(target_outcome$N, na.rm = TRUE)
  )
}

## =========================================================
## 5. Run one MRBEE-TL analysis
## =========================================================

wald_p <- function(beta, se) {
  2 * pnorm(abs(beta / se), lower.tail = FALSE)
}

run_mrbee_tl <- function(dat, auxiliary, target, exposure, outcome, threshold) {
  exp_dat <- dat$exposure
  C_source <- dat$background$C_source
  C_target <- dat$background$C_target

  # 1. Select candidate IVs using the joint FEMA exposure chi-square statistic.
  fema_chi2 <- exp_dat$z_source^2 / C_source[1, 1] +
    exp_dat$z_target^2 / C_target[1, 1]
  fema_p <- pchisq(fema_chi2, df = 2, lower.tail = FALSE)
  candidates <- data.table(SNP = exp_dat$SNP, P = fema_p)[P < threshold]
  if (nrow(candidates) == 0) stop("No FEMA-significant SNP")

  # 2. Apply 10,000 kb LD clumping with the EUR reference.
  if (nrow(candidates) == 1) {
    fema_iv <- candidates$SNP
  } else {
    clumped <- ieugwasr::ld_clump(
      data.frame(rsid = candidates$SNP, pval = candidates$P),
      clump_kb = clump_kb,
      clump_r2 = clump_r2,
      clump_p = threshold,
      bfile = eur_bfile,
      plink_bin = plink_bin
    )
    fema_iv <- clumped$rsid
  }
  if (length(fema_iv) == 0) stop("No IV remains after clumping")

  # 3. Remove MR-GxE outliers among source genome-wide significant IVs.
  iv_dat <- exp_dat[match(fema_iv, exp_dat$SNP)]
  source_gws <- which(iv_dat$p_source < gxe_p_threshold)
  gxe_removed <- character(0)

  if (length(source_gws) >= 3) {
    C_gxe <- diag(c(C_source[1, 1], C_target[1, 1]))
    fit_gxe <- MRBEE::MRBEE.IMRP.UV(
      by = iv_dat$b_target[source_gws],
      bx = iv_dat$b_source[source_gws],
      byse = iv_dat$se_target[source_gws],
      bxse = iv_dat$se_source[source_gws],
      Rxy = C_gxe,
      adjust.method = "Sidak"
    )
    delta <- as.numeric(fit_gxe$delta)
    gxe_removed <- iv_dat$SNP[source_gws][which(is.finite(delta) & delta != 0)]
  }

  final_iv <- setdiff(fema_iv, gxe_removed)
  if (length(final_iv) < 3) stop("Fewer than three IVs after MR-GxE filtering")

  # 4. Attach and harmonize source and target outcome associations after IV selection.
  mr_dat <- add_outcomes(
    exp_dat, dat$source_outcome, dat$target_outcome, final_iv
  )
  if (nrow(mr_dat) < 3) stop("Fewer than three harmonised IVs")

  # 5. Estimate the source causal effect from source exposure and outcome data.
  fit_source <- MRBEE::MRBEE.IMRP.UV(
    by = mr_dat$b_out_source,
    bx = mr_dat$b_source,
    byse = mr_dat$se_out_source,
    bxse = mr_dat$se_source,
    Rxy = C_source
  )

  # 6. Estimate the target causal effect using target-population data.
  fit_target <- MRBEE::MRBEE.IMRP.UV(
    by = mr_dat$b_out_target,
    bx = mr_dat$b_target,
    byse = mr_dat$se_out_target,
    bxse = mr_dat$se_target,
    Rxy = C_target
  )

  # 7. Fit the single-exposure transfer model with the source MRBEE estimate as prior.
  fit_tl <- MRBEEX::MRBEE_TL_UV(
    by = mr_dat$b_out_target,
    bX = matrix(mr_dat$b_target, ncol = 1),
    byse = mr_dat$se_out_target,
    bXse = matrix(mr_dat$se_target, ncol = 1),
    Rxy = C_target,
    theta.source = as.numeric(fit_source$theta)[1],
    theta.source.cov = matrix(as.numeric(fit_source$vartheta)[1], nrow = 1),
    LD = "identity",
    pip.thres = 0.5,
    ebic.delta = 0,
    ebic.gamma = 1,
    reliability.thres = 0.5,
    sampling.time = sampling_time
  )

  # 8. Save source, target and transfer-learning estimates.
  result <- data.table(
    method = c(
      "MRBEE_source_specific",
      "MRBEE_target_specific",
      "MRBEE_TL_UV"
    ),
    population = c(auxiliary, target, target),
    beta = c(
      as.numeric(fit_source$theta)[1],
      as.numeric(fit_target$theta)[1],
      as.numeric(fit_tl$theta)[1]
    ),
    se = c(
      sqrt(as.numeric(fit_source$vartheta)[1]),
      sqrt(as.numeric(fit_target$vartheta)[1]),
      as.numeric(fit_tl$theta.se)[1]
    )
  )
  if (any(!is.finite(result$beta)) || any(!is.finite(result$se)) || any(result$se <= 0)) {
    stop("MRBEE-TL returned an invalid estimate")
  }
  result[, `:=`(
    analysis_id = paste(auxiliary, target, exposure, outcome, threshold, sep = "|"),
    auxiliary = auxiliary,
    target = target,
    exposure = exposure,
    outcome = outcome,
    threshold = threshold,
    z = beta / se,
    p = mapply(wald_p, beta, se),
    CI_low = beta - 1.96 * se,
    CI_high = beta + 1.96 * se,
    nIV = nrow(mr_dat),
    n_fema_candidates = nrow(candidates),
    n_fema_iv = length(fema_iv),
    n_gxe_removed = length(gxe_removed),
    source_exposure_N = dat$source_exposure_N,
    source_outcome_N = dat$source_outcome_N,
    target_exposure_N = dat$target_exposure_N,
    target_outcome_N = dat$target_outcome_N,
    clump_kb = clump_kb,
    clump_r2 = clump_r2,
    gxe_p_threshold = gxe_p_threshold,
    sampling_strategy = c(NA_character_, NA_character_, "half_sample_subsampling"),
    sampling_time = c(NA_integer_, NA_integer_, sampling_time),
    background_file = dat$background$filename,
    status = "ok",
    error = NA_character_
  )]

  setcolorder(result, c(
    "analysis_id", "method", "population", "auxiliary", "target",
    "exposure", "outcome", "threshold", "beta", "se", "z", "p",
    "CI_low", "CI_high", "nIV", "n_fema_candidates", "n_fema_iv",
    "n_gxe_removed", "source_exposure_N", "source_outcome_N",
    "target_exposure_N", "target_outcome_N", "clump_kb", "clump_r2",
    "gxe_p_threshold", "sampling_strategy", "sampling_time",
    "background_file", "status", "error"
  ))
  result
}

## =========================================================
## 6. Run analyses and append results
## =========================================================

make_failure <- function(auxiliary, target, exposure, outcome, threshold, error) {
  data.table(
    analysis_id = paste(auxiliary, target, exposure, outcome, threshold, sep = "|"),
    method = "MRBEE_TL_UV", population = target,
    auxiliary, target, exposure, outcome, threshold,
    beta = NA_real_, se = NA_real_, z = NA_real_, p = NA_real_,
    CI_low = NA_real_, CI_high = NA_real_, nIV = NA_integer_,
    n_fema_candidates = NA_integer_, n_fema_iv = NA_integer_,
    n_gxe_removed = NA_integer_, source_exposure_N = NA_real_,
    source_outcome_N = NA_real_, target_exposure_N = NA_real_,
    target_outcome_N = NA_real_, clump_kb, clump_r2, gxe_p_threshold,
    sampling_strategy = "half_sample_subsampling", sampling_time,
    background_file = NA_character_, status = "failed", error = as.character(error)[1]
  )
}

append_result <- function(result) {
  has_result <- file.exists(result_file) && file.info(result_file)$size > 0
  fwrite(
    result, result_file, sep = "\t", quote = FALSE, na = "NA",
    append = has_result, col.names = !has_result
  )
}

completed <- character(0)
if (file.exists(result_file) && file.info(result_file)$size > 0) {
  completed <- unique(fread(result_file, select = "analysis_id")$analysis_id)
}

for (i in seq_len(nrow(combinations))) {
  auxiliary <- combinations$auxiliary[i]
  target <- combinations$target[i]

  for (exposure in exposures) {
    for (outcome in outcomes) {
      message("Preparing: ", auxiliary, " -> ", target, " | ", exposure, " -> ", outcome)

      preparation_error <- NULL
      prepared <- tryCatch(
        prepare_pair(auxiliary, target, exposure, outcome),
        error = function(e) {
          preparation_error <<- conditionMessage(e)
          NULL
        }
      )

      for (threshold in thresholds) {
        analysis_id <- paste(auxiliary, target, exposure, outcome, threshold, sep = "|")
        if (analysis_id %in% completed) next

        result <- if (is.null(prepared)) {
          make_failure(auxiliary, target, exposure, outcome, threshold, preparation_error)
        } else {
          tryCatch(
            run_mrbee_tl(prepared, auxiliary, target, exposure, outcome, threshold),
            error = function(e) make_failure(
              auxiliary, target, exposure, outcome, threshold, conditionMessage(e)
            )
          )
        }

        append_result(result)
        completed <- c(completed, analysis_id)
        message("  ", threshold, " appended [", result$status[1], "]")
      }

      rm(prepared)
      invisible(gc())
    }
  }
}

message("MRBEE-TL finished: ", result_file)
