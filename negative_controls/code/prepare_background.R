#!/usr/bin/env Rscript

# Estimate the target-population backgrounds used by MRBEE and MRBEE-TL.
# Stage 2 retains the intercepts estimated in stage 1, including cross-population
# intercepts. The XMR runner estimates its own zero-cross-intercept backgrounds.
.xmr_args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(.xmr_args)) stop("Run this entry point with Rscript.")
.xmr_root <- dirname(normalizePath(sub("^--file=", "", .xmr_args[[1L]])))
while (!file.exists(file.path(.xmr_root, ".xmr-reproduce-root"))) {
  parent <- dirname(.xmr_root)
  if (identical(parent, .xmr_root)) stop("Repository root not found.")
  .xmr_root <- parent
}
setwd(.xmr_root)
suppressPackageStartupMessages({library(XMR); library(XMAP); library(data.table)})
if (!"cross_intercept_zero" %in% names(formals(XMR::estimate_bgparas))) {
  stop("Install the XMR version specified in environment/README.md.")
}
set.seed(1)
nc_dir <- file.path(.xmr_root, "negative_controls")
formatted_dir <- file.path(nc_dir, "formatted_data")
output_dir <- file.path(nc_dir, "results", "background_MRBEE_TL_target")
info <- read.csv(file.path(nc_dir, "info.csv"), stringsAsFactors = FALSE)
if (!"phenocode" %in% names(info)) stop("info.csv must contain phenocode.")
outcomes <- c("1747_5", "1727")
exposures <- setdiff(unique(as.character(info$phenocode)), outcomes)
if (!length(exposures) || anyNA(exposures) || any(!nzchar(exposures))) {
  stop("info.csv contains missing exposure codes.")
}

load_scores <- function(target) {
  if (target == "AFR") {
    ld <- fread(file.path(nc_dir, "LDscore_eas_brit_afr.txt"))
    ambiguous <- which(ld$allele1 == XMR::comple(ld$allele2))
    if (length(ambiguous)) ld <- ld[-ambiguous, ]
    ld$EUR_AFR <- ld$AFR_EUR
    ld$EAS_AFR <- ld$AFR_EAS
    ld$EAS_EUR <- ld$EUR_EAS
    ld$EUR_EUR <- ld$EUR
    ld$EAS_EAS <- ld$EAS
    ld$AFR_AFR <- ld$AFR
    snp_col <- "rsid"
  } else {
    ld <- read.csv(file.path(nc_dir, "LDscore_SAS_EUR.csv"), sep = "\t")
    ld$EUR_CSA <- ld$CSA_EUR
    ld$EUR_EUR <- ld$EUR
    ld$CSA_CSA <- ld$CSA
    snp_col <- "SNP"
  }
  list(data = as.data.frame(ld), snp_col = snp_col)
}

estimate_pair <- function(target, exposure, outcome, scores) {
  paras <- XMR::estimate_bgparas(
    ldscore0 = scores$data, exposure_code = exposure, outcome_code = outcome,
    pop_ref = "EUR", pop_target = target, gwas_dir = formatted_dir,
    snp_col = scores$snp_col, z2_threshold = 30, cross_intercept_zero = FALSE
  )
  columns <- c("trait1", "trait2", "pop1", "pop2", "omega11", "omega12",
               "omega22", "c1", "c12", "c2")
  if (is.null(paras) || nrow(paras) != 6L || !all(columns %in% names(paras))) {
    stop("Incomplete background for ", target, ": ", exposure, " -> ", outcome)
  }
  paras <- paras[, columns]
  if (anyDuplicated(paras[, c("trait1", "trait2")]) ||
      any(!is.finite(as.matrix(paras[, columns[5:10]])))) {
    stop("Invalid background for ", target, ": ", exposure, " -> ", outcome)
  }
  filename <- file.path(output_dir, paste0(target, "_", outcome), paste0("bgparas_", exposure))
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(filename, ".tmp")
  write.table(paras, temporary, quote = FALSE, col.names = FALSE, row.names = FALSE)
  if (!file.rename(temporary, filename)) stop("Cannot save background: ", filename)
  message("Saved ", target, ": ", exposure, " -> ", outcome)
}

for (target in c("AFR", "CSA")) {
  scores <- load_scores(target)
  for (outcome in outcomes) for (exposure in exposures) {
    estimate_pair(target, exposure, outcome, scores)
  }
}
message("Background parameters written to ", output_dir)
