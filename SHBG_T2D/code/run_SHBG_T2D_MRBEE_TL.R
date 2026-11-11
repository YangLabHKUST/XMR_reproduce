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
# MRBEE-TL analysis: SHBG (30830) -> T2D (250.2)
#
# Source: EUR
# Targets: EAS, AFR, CSA
# IV thresholds: 5e-5 and 5e-8
# LD clumping: 10000 kb, r2 = 0.001
#
# Instrument workflow:
#   1. Intersect SHBG SNPs across EUR/EAS/AFR/CSA.
#   2. Harmonise all SHBG effects to the EUR allele orientation.
#   3. Calculate the MRBEE-TL four-population joint-chi-square FEMA P.
#   4. Run two clumping modes while keeping the exposure-only FEMA candidates fixed:
#      - EUR_clump: one EUR-panel clump shared by source and target candidates.
#      - target_specific_clump: target/TL FEMA sets use matched LD panels.
#   5. Estimate the EUR source effect from FEMA-selected, EUR-LD-clumped IVs
#      with matched EUR SHBG/T2D associations, without an extra EUR P filter.
#   6. Match target outcomes separately; do not intersect outcome availability
#      across populations. Actual source and target IV sets may therefore differ.
#   7. Remove target-specific MR-GxE outliers from the target/TL IV set.
#   8. Run MRBEE_TL_UV with half-sample subsampling.
############################################################

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ieugwasr)
  library(MRBEE)
  library(MRBEEX)
})

set.seed(1)

get_arg <- function(args, key, default) {
  hit <- grep(paste0("^", key, "="), args)
  if (length(hit)) return(sub(paste0("^", key, "="), "", args[hit[1]]))
  hit <- which(args == key)
  if (length(hit) && hit[1] < length(args)) return(args[hit[1] + 1L])
  default
}

threshold_label <- function(x) {
  if (isTRUE(all.equal(x, 5e-5))) "5e-05" else "5e-08"
}

complement <- function(x) chartr("ACGT", "TGCA", toupper(x))

args <- commandArgs(trailingOnly = TRUE)
formatted_dir <- get_arg(
  args, "--formatted_dir",
  file.path(.xmr_root, "SHBG_T2D/formatted_data")
)
reference_dir <- get_arg(
  args, "--reference_dir",
  file.path(.xmr_root, "negative_controls")
)
background_file <- get_arg(
  args, "--background_file",
  file.path(.xmr_root, "SHBG_T2D/results/SHBG-T2D-bgparas.RData")
)
output_dir <- get_arg(
  args, "--output_dir",
  file.path(.xmr_root, "SHBG_T2D/results/MRBEE_TL_FEMA_source_results_clump10000kb")
)

exposure_code <- "30830"
outcome_code <- "250.2"
pops <- c("EUR", "EAS", "AFR", "CSA")
targets <- c("EAS", "AFR", "CSA")
thresholds <- c(5e-5, 5e-8)
clump_modes <- c("EUR_clump", "target_specific_clump")
clump_kb <- 10000L
clump_r2 <- 0.001
gxe_p_threshold <- 5e-8
sampling_time <- 1000L

bfiles <- c(
  EUR = file.path(reference_dir, "1kg_pops/EUR_phase3_hg38_autosomes"),
  EAS = file.path(reference_dir, "1kg_pops/EAS_phase3_hg38_autosomes"),
  AFR = file.path(reference_dir, "1kg_pops/AFR_phase3_hg38_autosomes"),
  CSA = file.path(reference_dir, "1kg_pops/SAS_phase3_hg38_autosomes")
)
plink_bin <- get_arg(args, "--plink_bin", file.path(reference_dir, "plink_Linux"))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
result_file <- file.path(output_dir, "SHBG_T2D_MRBEE_TL_all_results.tsv")
if (file.exists(result_file)) file.remove(result_file)

read_gwas <- function(pop, trait) {
  filename <- file.path(formatted_dir, paste0(pop, "_", trait))
  if (!file.exists(filename)) stop("Missing GWAS file: ", filename)
  dat <- fread(filename)
  required <- c("SNP", "A1", "A2", "Z", "N", "P")
  missing <- setdiff(required, names(dat))
  if (length(missing)) stop(filename, " is missing: ", paste(missing, collapse = ", "))

  dat[, `:=`(
    SNP = as.character(SNP), A1 = toupper(as.character(A1)),
    A2 = toupper(as.character(A2)), Z = as.numeric(Z),
    N = as.numeric(N), P = as.numeric(P)
  )]
  # Inputs are QC-filtered; construct the model fields and standardized effects.
  dat[, `:=`(se = 1 / sqrt(N), beta = Z / sqrt(N))]
  as.data.frame(dat)
}

# Align one population to the EUR exposure orientation, including strand flips.
align_to_reference <- function(reference, other) {
  other <- other[match(reference$SNP, other$SNP), , drop = FALSE]
  same <- reference$A1 == other$A1 & reference$A2 == other$A2
  swap <- reference$A1 == other$A2 & reference$A2 == other$A1
  strand_same <- reference$A1 == complement(other$A1) &
    reference$A2 == complement(other$A2)
  strand_swap <- reference$A1 == complement(other$A2) &
    reference$A2 == complement(other$A1)
  keep <- same | swap | strand_same | strand_swap
  flip <- !same & (swap | strand_swap)

  other$beta[flip] <- -other$beta[flip]
  other$Z[flip] <- -other$Z[flip]
  other$A1 <- reference$A1
  other$A2 <- reference$A2
  list(data = other, keep = keep)
}

shbg <- setNames(lapply(pops, read_gwas, trait = exposure_code), pops)

common_shbg <- Reduce(intersect, lapply(shbg, `[[`, "SNP"))
if (!length(common_shbg)) stop("No SHBG SNP is shared by all four populations")

shbg <- lapply(shbg, function(x) x[match(common_shbg, x$SNP), , drop = FALSE])
reference <- shbg$EUR
allele_keep <- rep(TRUE, nrow(reference))
for (pop in targets) {
  aligned <- align_to_reference(reference, shbg[[pop]])
  shbg[[pop]] <- aligned$data
  allele_keep <- allele_keep & aligned$keep
}
shbg <- lapply(shbg, function(x) x[allele_keep, , drop = FALSE])
reference <- shbg$EUR
message("Four-population allele-aligned SHBG SNPs: ", nrow(reference))

# MRBEE-TL needs the LDSC intercept matrices both for its FEMA joint test and
# for the downstream bias-corrected estimating equations.
if (!file.exists(background_file)) stop("Missing background file: ", background_file,
  "; run SHBG_T2D/code/prepare_background.R first.")
load(background_file)
if (!exists("C")) stop("Background file must contain C")
C_list <- C
rm(C)

# Follow the supplied MRBEE-TL implementation: combine the four ancestry
# exposure Z statistics as a joint chi-square and scale each component by its
# LDSC estimation-error variance. The shared-IV analysis then uses EUR LD.
fema_chi2 <- Reduce(`+`, lapply(pops, function(pop) {
  shbg[[pop]]$Z^2 / C_list[[pop]][2, 2]
}))
fema <- data.frame(
  SNP = reference$SNP, A1 = reference$A1, A2 = reference$A2,
  chi2_fema = fema_chi2, df_fema = length(pops),
  p_fema = pmax(
    pchisq(fema_chi2, df = length(pops), lower.tail = FALSE),
    .Machine$double.xmin
  )
)
fwrite(fema, file.path(output_dir, "SHBG_four_population_FEMA.tsv"), sep = "\t")

# Align SHBG and T2D independently within each population. Outcome availability
# is not intersected across all four populations and does not define FEMA IVs.
mr_data <- setNames(vector("list", length(pops)), pops)
for (pop in pops) {
  outcome <- read_gwas(pop, outcome_code)
  snps <- intersect(shbg[[pop]]$SNP, outcome$SNP)
  exposure <- shbg[[pop]][match(snps, shbg[[pop]]$SNP), , drop = FALSE]
  outcome <- outcome[match(snps, outcome$SNP), , drop = FALSE]
  aligned_outcome <- align_to_reference(exposure, outcome)
  keep <- aligned_outcome$keep
  outcome <- aligned_outcome$data

  mr_data[[pop]] <- data.frame(
    SNP = exposure$SNP[keep],
    A1 = exposure$A1[keep],
    A2 = exposure$A2[keep],
    bx = exposure$beta[keep],
    by = outcome$beta[keep],
    bxse = exposure$se[keep],
    byse = outcome$se[keep],
    px = exposure$P[keep],
    py = outcome$P[keep]
  )
}
fema_mr <- fema

fit_mrbee <- function(pop, snps, population_data = mr_data[[pop]]) {
  dat <- population_data[match(snps, population_data$SNP), , drop = FALSE]
  if (nrow(dat) < 3) stop("Fewer than three IVs for ", pop)
  fit <- MRBEE::MRBEE.IMRP.UV(
    by = as.numeric(dat$by), bx = as.numeric(dat$bx),
    byse = as.numeric(dat$byse), bxse = as.numeric(dat$bxse),
    Rxy = C_list[[pop]][2:3, 2:3]
  )
  list(fit = fit, dat = dat)
}

panel_label <- function(pop) {
  if (pop == "CSA") "SAS target panel" else paste0(pop, " panel")
}

result_row <- function(threshold, clump_mode, pop, target_pop, method,
                       n_iv, beta, se, clump_panel, shared_iv,
                       n_gxe_removed = NA_integer_,
                       iv_selection = "MRBEE-TL four-population joint-chi-square FEMA P",
                       status = "ok", error = NA_character_) {
  data.table(
    method = method, clump_mode = clump_mode, population = pop,
    target_population = target_pop,
    source_population = "EUR", exposure = exposure_code, outcome = outcome_code,
    threshold = threshold, beta = beta, se = se,
    p = if (is.finite(beta) && is.finite(se) && se > 0) {
      2 * pnorm(abs(beta / se), lower.tail = FALSE)
    } else {
      NA_real_
    },
    nIV = n_iv, n_gxe_removed = n_gxe_removed,
    iv_selection = iv_selection,
    source_iv_selection = "Four-population joint-chi-square FEMA P; EUR LD clump",
    shared_IV_across_populations = FALSE,
    # FEMA flags refer to clumped candidates, before matching and MR-GxE filtering.
    FEMA_IV_shared_across_populations = shared_iv,
    FEMA_IV_shared_across_targets = shared_iv,
    clump_kb = clump_kb, clump_r2 = clump_r2,
    clump_panel = clump_panel,
    gxe_p_threshold = gxe_p_threshold,
    sampling_strategy = if (method == "MRBEE_TL_UV") {
      "half_sample_subsampling"
    } else {
      NA_character_
    },
    sampling_time = if (method == "MRBEE_TL_UV") sampling_time else NA_integer_,
    status = status, error = error
  )
}

clump_fema <- function(candidates, pop, threshold) {
  clumped <- ieugwasr::ld_clump(
    data.frame(rsid = candidates$SNP, pval = candidates$p_fema),
    clump_kb = clump_kb, clump_r2 = clump_r2, clump_p = threshold,
    bfile = bfiles[[pop]], plink_bin = plink_bin
  )
  if (is.null(clumped) || !nrow(clumped)) return(character(0))
  intersect(as.character(clumped$rsid), fema_mr$SNP)
}

run_mode <- function(threshold, clump_mode) {
  label <- threshold_label(threshold)
  shared_iv <- clump_mode == "EUR_clump"
  candidates <- fema_mr[fema_mr$p_fema < threshold, , drop = FALSE]
  if (!nrow(candidates)) {
    return(result_row(
      threshold, clump_mode, "ALL", "ALL", "MRBEE_TL_UV", 0L,
      NA_real_, NA_real_, "not_run", shared_iv,
      status = "failed", error = "No FEMA-significant SNP"
    ))
  }

  if (shared_iv) {
    eur_ivs <- clump_fema(candidates, "EUR", threshold)
    ivs_by_pop <- setNames(rep(list(eur_ivs), length(pops)), pops)
  } else {
    ivs_by_pop <- setNames(
      lapply(pops, function(pop) clump_fema(candidates, pop, threshold)),
      pops
    )
  }

  # Reuse the EUR FEMA clump, matching the source outcome only after selection.
  # Target outcome availability and MR-GxE exclusions do not restrict this fit.
  source_ivs <- intersect(ivs_by_pop$EUR, mr_data$EUR$SNP)
  source_dat <- mr_data$EUR[match(source_ivs, mr_data$EUR$SNP), , drop = FALSE]
  source_dat$p_fema <- fema_mr$p_fema[match(source_ivs, fema_mr$SNP)]
  fwrite(
    source_dat,
    file.path(output_dir, paste0("MRBEE_source_EUR_FEMA_", label, "_IVs.tsv")),
    sep = "\t", quote = FALSE
  )
  source_fit <- tryCatch(
    fit_mrbee("EUR", source_ivs),
    error = function(e) e
  )

  # Save clumped candidates before outcome matching and target MR-GxE filtering.
  for (pop in pops) {
    ivs <- ivs_by_pop[[pop]]
    fwrite(
      fema_mr[fema_mr$SNP %in% ivs, ],
      file.path(
        output_dir,
        paste0("MRBEE_TL_", clump_mode, "_", pop, "_", label, "_IVs.tsv")
      ),
      sep = "\t", quote = FALSE
    )
  }

  rows <- list()
  j <- 1L
  for (pop in targets) {
    clumped_ivs <- ivs_by_pop[[pop]]
    target_panel <- if (shared_iv) "EUR panel" else panel_label(pop)

    # Record the same FEMA-based EUR source estimate for each target comparison.
    if (inherits(source_fit, "error")) {
      rows[[j]] <- result_row(
        threshold, clump_mode, "EUR", pop, "MRBEE_source_specific",
        length(source_ivs), NA_real_, NA_real_, "EUR panel", shared_iv,
        n_gxe_removed = 0L,
        status = "failed", error = conditionMessage(source_fit)
      )
    } else {
      theta_source <- as.numeric(source_fit$fit$theta)[1]
      var_source <- as.numeric(source_fit$fit$vartheta)[1]
      rows[[j]] <- result_row(
        threshold, clump_mode, "EUR", pop, "MRBEE_source_specific",
        length(source_ivs), theta_source, sqrt(var_source), "EUR panel", shared_iv,
        n_gxe_removed = 0L
      )
    }
    j <- j + 1L

    # Target/TL FEMA IVs need target outcome and EUR exposure for MR-GxE, but
    # do not require EUR outcome availability or membership in source_ivs.
    clumped_ivs <- Reduce(intersect, list(
      clumped_ivs, shbg$EUR$SNP, mr_data[[pop]]$SNP
    ))
    source_exp <- shbg$EUR[match(clumped_ivs, shbg$EUR$SNP), , drop = FALSE]
    target_exp <- mr_data[[pop]][match(clumped_ivs, mr_data[[pop]]$SNP), , drop = FALSE]

    # Match the other two workflows: detect MR-GxE outliers among EUR
    # genome-wide significant IVs by treating target exposure as the response.
    source_gws <- which(source_exp$P < gxe_p_threshold)
    gxe_removed <- character(0)
    if (length(source_gws) >= 3L) {
      C_gxe <- diag(c(C_list$EUR[2, 2], C_list[[pop]][2, 2]))
      fit_gxe <- MRBEE::MRBEE.IMRP.UV(
        by = target_exp$bx[source_gws],
        bx = source_exp$beta[source_gws],
        byse = target_exp$bxse[source_gws],
        bxse = source_exp$se[source_gws],
        Rxy = C_gxe,
        adjust.method = "Sidak"
      )
      delta <- as.numeric(fit_gxe$delta)
      gxe_removed <- source_exp$SNP[source_gws][
        which(is.finite(delta) & delta != 0)
      ]
    }

    final_ivs <- setdiff(clumped_ivs, gxe_removed)
    if (length(final_ivs) < 3L) {
      error <- "Fewer than three IVs after MR-GxE filtering"
      for (spec in list(
        c(pop, "MRBEE_target_specific"),
        c(pop, "MRBEE_TL_UV")
      )) {
        rows[[j]] <- result_row(
          threshold, clump_mode, spec[1], pop, spec[2], length(final_ivs),
          NA_real_, NA_real_, target_panel, shared_iv,
          n_gxe_removed = length(gxe_removed), status = "failed", error = error
        )
        j <- j + 1L
      }
      next
    }

    fwrite(
      fema_mr[fema_mr$SNP %in% final_ivs, ],
      file.path(
        output_dir,
        paste0(
          "MRBEE_TL_", clump_mode, "_", pop, "_", label,
          "_post_MR_GxE_IVs.tsv"
        )
      ),
      sep = "\t", quote = FALSE
    )

    target_fit <- tryCatch(fit_mrbee(pop, final_ivs), error = function(e) e)

    if (inherits(target_fit, "error")) {
      error <- conditionMessage(target_fit)
      for (spec in list(
        c(pop, "MRBEE_target_specific"),
        c(pop, "MRBEE_TL_UV")
      )) {
        rows[[j]] <- result_row(
          threshold, clump_mode, spec[1], pop, spec[2], length(final_ivs),
          NA_real_, NA_real_, target_panel, shared_iv,
          n_gxe_removed = length(gxe_removed), status = "failed", error = error
        )
        j <- j + 1L
      }
      next
    }

    specific_beta <- as.numeric(target_fit$fit$theta)[1]
    specific_se <- sqrt(as.numeric(target_fit$fit$vartheta)[1])
    rows[[j]] <- result_row(
      threshold, clump_mode, pop, pop, "MRBEE_target_specific",
      length(final_ivs), specific_beta, specific_se, target_panel,
      shared_iv, n_gxe_removed = length(gxe_removed)
    )
    j <- j + 1L

    if (inherits(source_fit, "error")) {
      rows[[j]] <- result_row(
        threshold, clump_mode, pop, pop, "MRBEE_TL_UV", length(final_ivs),
        NA_real_, NA_real_, target_panel, shared_iv,
        n_gxe_removed = length(gxe_removed), status = "failed",
        error = paste("FEMA-based EUR source fit failed:", conditionMessage(source_fit))
      )
      j <- j + 1L
      next
    }

    dat <- target_fit$dat
    tl_fit <- tryCatch(
      MRBEEX::MRBEE_TL_UV(
        by = as.numeric(dat$by),
        bX = matrix(as.numeric(dat$bx), ncol = 1),
        byse = as.numeric(dat$byse),
        bXse = matrix(as.numeric(dat$bxse), ncol = 1),
        Rxy = C_list[[pop]][2:3, 2:3],
        theta.source = theta_source,
        theta.source.cov = matrix(var_source, 1, 1),
        LD = "identity",
        pip.thres = 0.5,
        ebic.delta = 0,
        ebic.gamma = 1,
        reliability.thres = 0.5,
        sampling.time = sampling_time
      ),
      error = function(e) e
    )

    if (inherits(tl_fit, "error")) {
      rows[[j]] <- result_row(
        threshold, clump_mode, pop, pop, "MRBEE_TL_UV", length(final_ivs),
        NA_real_, NA_real_, target_panel, shared_iv,
        n_gxe_removed = length(gxe_removed),
        status = "failed", error = conditionMessage(tl_fit)
      )
    } else {
      rows[[j]] <- result_row(
        threshold, clump_mode, pop, pop, "MRBEE_TL_UV", length(final_ivs),
        as.numeric(tl_fit$theta)[1], as.numeric(tl_fit$theta.se)[1],
        target_panel, shared_iv, n_gxe_removed = length(gxe_removed)
      )
    }
    j <- j + 1L
  }

  mode_results <- rbindlist(rows, fill = TRUE)
  mode_results
}

all_results <- list()
k <- 1L
for (threshold in thresholds) {
  for (clump_mode in clump_modes) {
    message(
      "MRBEE-TL | ", clump_mode,
      " | threshold ", threshold_label(threshold)
    )
    all_results[[k]] <- run_mode(threshold, clump_mode)
    k <- k + 1L
  }
}

final <- rbindlist(all_results, fill = TRUE)
setorder(final, threshold, clump_mode, population, method)
fwrite(final, result_file, sep = "\t", quote = FALSE, na = "NA")
message("MRBEE-TL finished: ", result_file)
