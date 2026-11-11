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
# Negative-control MRAPSS: EUR-only analysis
#
# Analysis design:
#   exposure: all phenocodes in info.csv except 1747_5 and 1727
#   outcome : 1747_5 and 1727
#   IV/MR   : EUR exposure/outcome data
#   IV thresholds: 5e-5 and 5e-8.
#   clumping ：10,000 kb，r2 = 0.001，EUR reference panel
#   Background: two-stage EUR LDSC for each exposure-outcome pair.
#
# Inputs have undergone MAF, INFO, duplicate-SNP and
# palindromic-SNP QC. Check model fields, align alleles,
# and validate finite values for LDSC and MR.
############################################################

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(XMR)
  library(XMAP)
  library(MRAPSS)
})

set.seed(1)

# 1. Set input paths, traits and analysis parameters.
nc_dir <- file.path(.xmr_root, "negative_controls")
formatted_dir <- file.path(nc_dir, "formatted_data")
info_file <- file.path(nc_dir, "info.csv")
result_dir <- file.path(.xmr_root, "negative_controls/results")
background_dir <- file.path(result_dir, "MRAPSS_EUR_background")
ldscore_file <- file.path(nc_dir, "LDscore_eas_brit_afr.txt")

plink_bin <- file.path(nc_dir, "plink_Linux")
bfiles <- c(
  EUR = file.path(nc_dir, "1kg_pops/EUR_phase3_hg38_autosomes")
)

targets <- "EUR"
thresholds <- c(5e-5, 5e-8)
clump_kb <- 10000
clump_r2 <- 0.001

result_file <- file.path(result_dir, "NC_MRAPSS_EUR_clump10000kb.tsv")

# info.csv defines the negative-control trait pairs.
if (!file.exists(info_file)) stop("Negative-control info file not found: ", info_file)
info <- fread(info_file)
if (!"phenocode" %in% names(info)) stop("info.csv must contain phenocode")
outcomes <- c("1747_5", "1727")
exposures <- setdiff(unique(as.character(info$phenocode)), outcomes)
exposures <- exposures[!is.na(exposures) & nzchar(exposures)]
if (length(exposures) == 0L) stop("No negative-control exposures found")

# 2. Check input directories and reference files; create the output directory.
if (!dir.exists(formatted_dir)) stop("Formatted-data directory not found: ", formatted_dir)
if (!file.exists(plink_bin)) stop("PLINK binary not found: ", plink_bin)

required_bfiles <- unlist(lapply(bfiles, paste0, c(".bed", ".bim", ".fam")))
missing_bfiles <- required_bfiles[!file.exists(required_bfiles)]
if (length(missing_bfiles) > 0) {
  stop("Missing target-population reference files: ", paste(missing_bfiles, collapse = ", "))
}

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(background_dir, recursive = TRUE, showWarnings = FALSE)

# 3. Define shared helper functions.
threshold_label <- function(x) {
  if (isTRUE(all.equal(x, 5e-5))) return("5e-05")
  if (isTRUE(all.equal(x, 5e-8))) return("5e-08")
  format(x, scientific = TRUE)
}

complement_allele <- function(x) chartr("ATCG", "TAGC", toupper(x))

scalar_or_na <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_real_)
  value <- suppressWarnings(as.numeric(x[1]))
  if (length(value) == 0 || !is.finite(value)) NA_real_ else value
}

# 4. Read extensionless GWAS files in the MR-APSS formatted schema.
read_gwas <- function(dataset, trait) {
  filename <- file.path(formatted_dir, paste0(dataset, "_", trait))
  if (!file.exists(filename)) stop("GWAS file not found: ", filename)

  dat <- fread(filename)
  required <- c("SNP", "A1", "A2", "Z", "N", "P")
  missing_cols <- setdiff(required, names(dat))
  if (length(missing_cols) > 0) {
    stop(filename, " is missing columns: ", paste(missing_cols, collapse = ", "))
  }

  # Inputs are QC-filtered; normalize field types while retaining the supplied rows.
  dat[, `:=`(
    SNP = as.character(SNP),
    A1 = toupper(as.character(A1)),
    A2 = toupper(as.character(A2)),
    Z = suppressWarnings(as.numeric(Z)),
    N = suppressWarnings(as.numeric(N)),
    P = suppressWarnings(as.numeric(P))
  )]
  dat[, .(SNP, A1, A2, Z, N, P)]
}

# 5. Read EUR LD scores for two-stage LDSC and the SNP/L2 input to MR-APSS.
load_ldscore <- function() {
  ld <- fread(ldscore_file)
  required <- c("rsid", "allele1", "allele2", "EUR")
  missing <- setdiff(required, names(ld))
  if (length(missing)) {
    stop(ldscore_file, " is missing: ", paste(missing, collapse = ", "))
  }
  ld[, `:=`(
    rsid = as.character(rsid),
    allele1 = toupper(as.character(allele1)),
    allele2 = toupper(as.character(allele2)),
    EUR = as.numeric(EUR)
  )]
  ld <- ld[!duplicated(rsid) & allele1 != complement_allele(allele2)]
  ld[, EUR_EUR := EUR]
  ld
}

ldscore_full <- load_ldscore()
ldscore_mrapss <- ldscore_full[is.finite(EUR), .(SNP = rsid, L2 = EUR)]

# Align GWAS Z scores to allele1 in the LD-score table.
align_z_to_ldscore <- function(dat, ld) {
  same <-
    (dat$A1 == ld$allele1 & dat$A2 == ld$allele2) |
    (dat$A1 == complement_allele(ld$allele1) &
       dat$A2 == complement_allele(ld$allele2))
  reverse <-
    (dat$A1 == ld$allele2 & dat$A2 == ld$allele1) |
    (dat$A1 == complement_allele(ld$allele2) &
       dat$A2 == complement_allele(ld$allele1))
  same[is.na(same)] <- FALSE
  reverse[is.na(reverse)] <- FALSE
  z <- dat$Z
  z[!(same | reverse)] <- NA_real_
  z[reverse] <- -z[reverse]
  z
}

# Cache each trait diagonal within this run for reuse across trait pairs.
ldsc_cache <- new.env(parent = emptyenv())

estimate_ldsc_entry <- function(dataset, trait1, trait2) {
  cache_key <- paste(dataset, trait1, trait2, sep = "|")
  reverse_key <- paste(dataset, trait2, trait1, sep = "|")
  if (exists(cache_key, ldsc_cache, inherits = FALSE)) {
    return(get(cache_key, ldsc_cache, inherits = FALSE))
  }
  if (exists(reverse_key, ldsc_cache, inherits = FALSE)) {
    return(get(reverse_key, ldsc_cache, inherits = FALSE))
  }

  dat1 <- read_gwas(dataset, trait1)
  dat2 <- read_gwas(dataset, trait2)
  common <- Reduce(intersect, list(ldscore_full$rsid, dat1$SNP, dat2$SNP))
  if (length(common) < 100L) stop("Too few EUR LDSC SNPs: ", trait1, " / ", trait2)

  dat1 <- dat1[match(common, SNP)]
  dat2 <- dat2[match(common, SNP)]
  ld <- ldscore_full[match(common, rsid)]
  z1 <- align_z_to_ldscore(dat1, ld)
  z2 <- align_z_to_ldscore(dat2, ld)

  usable <-
    is.finite(z1) & is.finite(z2) &
    is.finite(dat1$N) & dat1$N > 0 &
    is.finite(dat2$N) & dat2$N > 0 & is.finite(ld$EUR)
  stage1 <- usable & z1^2 < 30 & z2^2 < 30
  if (sum(stage1) < 100L || sum(usable) < 100L) {
    stop("Too few usable EUR LDSC SNPs: ", trait1, " / ", trait2)
  }
  regression_weight <- 1 / pmax(ld$EUR, 1)

  # Stage 1: estimate LDSC intercepts using SNPs below the Z-squared threshold.
  fit_stage1 <- XMAP::estimate_gc(
    data.frame(Z = z1[stage1], N = dat1$N[stage1]),
    data.frame(Z = z2[stage1], N = dat2$N[stage1]),
    ld$EUR[stage1], ld$EUR[stage1], ld$EUR_EUR[stage1],
    reg_w1 = regression_weight[stage1],
    reg_w2 = regression_weight[stage1],
    reg_wx = regression_weight[stage1],
    constrain_intercept = FALSE
  )
  intercepts <- c(
    as.numeric(fit_stage1$tau1$coefs[1]),
    as.numeric(fit_stage1$tau2$coefs[1]),
    as.numeric(fit_stage1$theta$coefs[1])
  )
  if (any(!is.finite(intercepts))) stop("Non-finite first-stage EUR LDSC intercept")

  # Stage 2: fix the intercepts and estimate Omega using all valid SNPs.
  fit_stage2 <- XMAP::estimate_gc(
    data.frame(Z = z1[usable], N = dat1$N[usable]),
    data.frame(Z = z2[usable], N = dat2$N[usable]),
    ld$EUR[usable], ld$EUR[usable], ld$EUR_EUR[usable],
    reg_w1 = regression_weight[usable],
    reg_w2 = regression_weight[usable],
    reg_wx = regression_weight[usable],
    constrain_intercept = TRUE,
    intercepts[1], intercepts[2], intercepts[3]
  )
  estimates <- c(
    Omega11 = as.numeric(fit_stage2$tau1$coefs[2]),
    Omega12 = as.numeric(fit_stage2$theta$coefs[2]),
    Omega22 = as.numeric(fit_stage2$tau2$coefs[2]),
    C11 = as.numeric(fit_stage2$tau1$coefs[1]),
    C12 = as.numeric(fit_stage2$theta$coefs[1]),
    C22 = as.numeric(fit_stage2$tau2$coefs[1])
  )
  if (any(!is.finite(estimates))) stop("Non-finite second-stage EUR LDSC estimate")
  assign(cache_key, estimates, envir = ldsc_cache)
  estimates
}

# Estimate the 2 x 2 C and Omega matrices for each EUR exposure-outcome pair.
get_background <- function(dataset, exposure, outcome) {
  message("    estimate EUR background: ", exposure, " -> ", outcome)
  xx <- estimate_ldsc_entry(dataset, exposure, exposure)
  xy <- estimate_ldsc_entry(dataset, exposure, outcome)
  yy <- estimate_ldsc_entry(dataset, outcome, outcome)
  traits <- c(paste0(dataset, "_", exposure), paste0(dataset, "_", outcome))
  Omega <- matrix(
    c(xx[["Omega11"]], xy[["Omega12"]], xy[["Omega12"]], yy[["Omega11"]]),
    2L, 2L, byrow = TRUE, dimnames = list(traits, traits)
  )
  C <- matrix(
    c(xx[["C11"]], xy[["C12"]], xy[["C12"]], yy[["C11"]]),
    2L, 2L, byrow = TRUE, dimnames = list(traits, traits)
  )
  if (any(!is.finite(C)) || any(diag(C) <= 0)) stop("Invalid EUR LDSC C matrix")
  if (any(!is.finite(Omega)) || any(diag(Omega) <= 0)) stop("Invalid EUR LDSC Omega matrix")

  rho <- Omega[1, 2] / sqrt(Omega[1, 1] * Omega[2, 2])
  if (is.finite(rho) && abs(rho) > 1) {
    Omega[1, 2] <- Omega[2, 1] <-
      0.95 * sign(rho) * sqrt(Omega[1, 1] * Omega[2, 2])
  }
  tryCatch(chol(C), error = function(e) stop("Estimated EUR C matrix is not positive definite"))

  filename <- file.path(background_dir, paste0(dataset, "_", exposure, "_", outcome, ".rds"))
  background <- list(C = C, Omega = Omega, file = filename)
  saveRDS(background, filename)
  background
}

# 9. Align exposure and outcome within the target population and attach L2.
prepare_pair <- function(dataset, exposure, outcome) {
  exposure_gwas <- read_gwas(dataset, exposure)
  outcome_gwas <- read_gwas(dataset, outcome)

  aligned <- XMR::align_data(
    dat1 = as.data.frame(exposure_gwas),
    dat2 = as.data.frame(outcome_gwas),
    trait1.name = paste0(dataset, "_", exposure),
    trait2.name = paste0(dataset, "_", outcome),
    h2.fix.intercept = FALSE,
    LDSC = FALSE
  )$dat

  if (ncol(aligned) != 9) stop("Unexpected XMR::align_data output")
  aligned <- as.data.table(aligned)
  setnames(aligned, c(
    "SNP", "A1", "A2", "b.exp", "b.out", "se.exp",
    "se.out", "pval.exp", "pval.out"
  ))

  aligned <- merge(aligned, ldscore_mrapss, by = "SNP")
  if (nrow(aligned) == 0) stop("No harmonised SNP overlaps the EUR LD-score file")

  list(
    data = aligned,
    background = get_background(dataset, exposure, outcome),
    exposure_N = median(exposure_gwas$N, na.rm = TRUE),
    outcome_N = median(outcome_gwas$N, na.rm = TRUE)
  )
}

# 10. Format each successful or failed analysis as one result row.
make_result <- function(dataset, exposure, outcome, threshold, status,
                        error = NA_character_, fit = NULL, prepared = NULL,
                        n_candidates = NA_integer_, n_iv = NA_integer_) {
  bg <- if (is.null(prepared)) NULL else prepared$background

  data.table(
    analysis_id = paste(dataset, exposure, outcome, threshold_label(threshold), sep = "|"),
    method = "MRAPSS",
    target_dataset = dataset,
    exposure = exposure,
    outcome = outcome,
    threshold = threshold,
    beta = if (is.null(fit)) NA_real_ else scalar_or_na(fit$beta),
    beta_se = if (is.null(fit)) NA_real_ else scalar_or_na(fit$beta.se),
    beta_pvalue = if (is.null(fit)) NA_real_ else scalar_or_na(fit$pvalue),
    tau_sq = if (is.null(fit)) NA_real_ else scalar_or_na(fit$tau.sq),
    pi0 = if (is.null(fit)) NA_real_ else scalar_or_na(fit$pi0),
    IVsignal_sum = if (is.null(fit)) NA_real_ else scalar_or_na(fit$post$IVsignal.sum),
    n_candidates = n_candidates,
    n_IV = n_iv,
    exposure_N = if (is.null(prepared)) NA_real_ else prepared$exposure_N,
    outcome_N = if (is.null(prepared)) NA_real_ else prepared$outcome_N,
    C11 = if (is.null(bg)) NA_real_ else bg$C[1, 1],
    C12 = if (is.null(bg)) NA_real_ else bg$C[1, 2],
    C22 = if (is.null(bg)) NA_real_ else bg$C[2, 2],
    Omega11 = if (is.null(bg)) NA_real_ else bg$Omega[1, 1],
    Omega12 = if (is.null(bg)) NA_real_ else bg$Omega[1, 2],
    Omega22 = if (is.null(bg)) NA_real_ else bg$Omega[2, 2],
    Cor_SelectionBias = TRUE,
    clump_kb = clump_kb,
    clump_r2 = clump_r2,
    background_file = if (is.null(bg)) NA_character_ else bg$file,
    status = status,
    error = as.character(error)[1]
  )
}

# 11. Select target-exposure IVs, clump and fit MR-APSS at one threshold.
run_mrapss <- function(prepared, dataset, exposure, outcome, threshold) {
  dat <- prepared$data
  n_candidates <- sum(is.finite(dat$pval.exp) & dat$pval.exp < threshold)
  if (n_candidates == 0) stop("No target-exposure IV candidate")

  clumped <- MRAPSS::clump(
    as.data.frame(dat),
    IV.Threshold = threshold,
    SNP_col = "SNP",
    pval_col = "pval.exp",
    clump_kb = clump_kb,
    clump_r2 = clump_r2,
    bfile = bfiles[[dataset]],
    plink_bin = plink_bin
  )

  n_iv <- if (is.null(clumped)) 0L else nrow(clumped)
  if (n_iv < 5L) stop("MRAPSS requires at least five IVs; found ", n_iv)

  clumped$Threshold <- threshold
  mr_data <- as.data.frame(clumped)[, c(
    "SNP", "A1", "A2", "b.exp", "se.exp", "pval.exp",
    "b.out", "se.out", "pval.out", "L2", "Threshold"
  )]

  fit <- MRAPSS::MRAPSS(
    mr_data,
    C = prepared$background$C,
    Omega = prepared$background$Omega,
    Cor.SelectionBias = TRUE
  )

  if (!is.finite(scalar_or_na(fit$beta)) || !is.finite(scalar_or_na(fit$beta.se))) {
    stop("MRAPSS returned a non-finite estimate")
  }

  make_result(
    dataset, exposure, outcome, threshold,
    status = "ok", fit = fit, prepared = prepared,
    n_candidates = n_candidates, n_iv = n_iv
  )
}

# 12. Read completed analysis IDs to support resuming the run.
results_started <- file.exists(result_file) && file.info(result_file)$size > 0
completed <- character(0)
if (results_started) {
  old_results <- fread(result_file)
  if (!"analysis_id" %in% names(old_results)) stop("Existing result file has no analysis_id")
  completed <- unique(old_results$analysis_id)
}

# 13. Append each completed threshold analysis to the result table.
append_result <- function(row) {
  fwrite(
    row, result_file, sep = "\t", quote = FALSE, na = "NA",
    append = results_started, col.names = !results_started
  )
  results_started <<- TRUE
}

# 14. Run all exposure-outcome pairs at both thresholds.
message(
  "MRAPSS analyses: ",
  length(targets) * length(exposures) * length(outcomes) * length(thresholds)
)
message("Result file: ", result_file)

for (dataset in targets) {
  for (exposure in exposures) {
    for (outcome in outcomes) {
      message("[", dataset, "] ", exposure, " -> ", outcome)

      preparation_error <- NA_character_
      prepared <- tryCatch(
        prepare_pair(dataset, exposure, outcome),
        error = function(e) {
          preparation_error <<- conditionMessage(e)
          NULL
        }
      )

      for (threshold in thresholds) {
        analysis_id <- paste(dataset, exposure, outcome, threshold_label(threshold), sep = "|")
        if (analysis_id %in% completed) {
          message("  skip existing ", threshold_label(threshold))
          next
        }

        result <- if (is.null(prepared)) {
          make_result(
            dataset, exposure, outcome, threshold,
            status = "failed", error = preparation_error
          )
        } else {
          tryCatch(
            run_mrapss(prepared, dataset, exposure, outcome, threshold),
            error = function(e) make_result(
              dataset, exposure, outcome, threshold,
              status = "failed", error = conditionMessage(e), prepared = prepared,
              n_candidates = sum(
                is.finite(prepared$data$pval.exp) &
                  prepared$data$pval.exp < threshold
              )
            )
          )
        }

        append_result(result)
        completed <- c(completed, analysis_id)
        message("  appended ", threshold_label(threshold), " [", result$status, "]")
      }

      rm(prepared)
      invisible(gc(verbose = FALSE))
    }
  }
}

message("Done: ", result_file)
