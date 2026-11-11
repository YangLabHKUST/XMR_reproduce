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
# Negative-control CAMeRa: EUR auxiliary population -> AFR/CSA targets
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
  library(GenomicRanges)
  library(IRanges)
  library(ieugwasr)
  library(CAMeRa)
})

set.seed(1)

# CAMeRa's regional FEMA helper uses sapply() to assemble SNP-by-population
# matrices. For a region containing exactly one common SNP, sapply() drops the
# matrix dimension and rowSums() fails. Preserve the official calculation while
# forcing that edge case back to a 1 x number-of-populations matrix.
fixed_effects_meta_analysis_fast_safe <- function(beta_mat, se_mat) {
  if (is.null(dim(beta_mat))) {
    beta_mat <- matrix(beta_mat, nrow = 1L)
  }
  if (is.null(dim(se_mat))) {
    se_mat <- matrix(se_mat, nrow = 1L)
  }

  weights <- 1 / se_mat^2
  beta <- rowSums(beta_mat * weights) / rowSums(weights, na.rm = TRUE)
  se <- sqrt(1 / rowSums(weights, na.rm = TRUE))
  z <- abs(beta / se)
  p <- pnorm(z, lower.tail = FALSE)
  nstudy <- apply(beta_mat, 1, function(x) sum(!is.na(x)))

  tibble::tibble(nstudy = nstudy, p = p, z = z)
}

assignInNamespace(
  "fixed_effects_meta_analysis_fast",
  fixed_effects_meta_analysis_fast_safe,
  ns = "CAMeRa"
)

## =========================================================
## 1. Analysis parameters
## =========================================================

nc_dir <- file.path(.xmr_root, "negative_controls")
formatted_dir <- file.path(nc_dir, "formatted_data")
info_file <- file.path(nc_dir, "info.csv")
result_dir <- file.path(.xmr_root, "negative_controls/results")

# Reference paths follow the input layout in DATA.md.
plink_bin <- file.path(nc_dir, "plink_Linux")
bfiles <- c(
  EUR = file.path(nc_dir, "1kg_pops/EUR_phase3_hg38_autosomes"),
  AFR = file.path(nc_dir, "1kg_pops/AFR_phase3_hg38_autosomes"),
  CSA = file.path(nc_dir, "1kg_pops/SAS_phase3_hg38_autosomes")
)

combinations <- data.table(
  auxiliary = c("EUR", "EUR"),
  target = c("AFR", "CSA")
)

thresholds <- c(5e-5, 5e-8)
clump_kb <- 10000
clump_r2 <- 0.001
region_radius <- 250000

# FALSE resumes existing results; TRUE starts a new result table.
overwrite_results <- FALSE
result_file <- file.path(result_dir, "NC_CAMeRa_EUR_to_AFR_CSA_clump10000kb.tsv")

# Read the negative-control trait pairs from info.csv.
if (!file.exists(info_file)) stop("Negative-control info file not found: ", info_file)
info <- fread(info_file)
if (!"phenocode" %in% names(info)) stop("info.csv must contain phenocode")
outcomes <- c("1747_5", "1727")
exposures <- setdiff(unique(as.character(info$phenocode)), outcomes)
exposures <- exposures[!is.na(exposures) & nzchar(exposures)]
if (length(exposures) == 0L) stop("No negative-control exposures found")

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
if (overwrite_results && file.exists(result_file)) file.remove(result_file)
if (!file.exists(plink_bin)) stop("PLINK binary not found: ", plink_bin)
required_bfiles <- unlist(lapply(bfiles, paste0, c(".bed", ".bim", ".fam")))
missing_bfiles <- required_bfiles[!file.exists(required_bfiles)]
if (length(missing_bfiles) > 0L) {
  stop("Missing LD reference files: ", paste(missing_bfiles, collapse = ", "))
}

## =========================================================
## 2. Read and harmonize GWAS data
## =========================================================

# Select the LD reference panel for each population.
ancestry_of <- function(prefix) {
  if (!prefix %in% names(bfiles)) stop("Unsupported population: ", prefix)
  prefix
}

# Formatted files may be extensionless or use .txt, .tsv or .gz.
find_gwas_file <- function(prefix, trait) {
  stem <- file.path(formatted_dir, paste0(prefix, "_", trait))
  files <- c(stem, paste0(stem, ".txt"), paste0(stem, ".tsv"), paste0(stem, ".gz"))
  files <- files[file.exists(files)]
  if (length(files) == 0) stop("Missing GWAS file: ", stem)
  files[1]
}

# Load each population BIM file once.
bim_cache <- new.env(parent = emptyenv())

read_bim <- function(prefix) {
  ancestry <- ancestry_of(prefix)
  if (exists(ancestry, bim_cache, inherits = FALSE)) {
    return(get(ancestry, bim_cache, inherits = FALSE))
  }

  bim <- fread(paste0(bfiles[[ancestry]], ".bim"), header = FALSE)
  setnames(bim, names(bim)[1:6], c("chr", "SNP", "cm", "position", "ref1", "ref2"))
  bim <- unique(bim[, .(SNP = as.character(SNP), chr, position)], by = "SNP")
  assign(ancestry, bim, envir = bim_cache)
  bim
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
  # Inputs are QC-filtered; convert the fields required by the model.
  dat[, `:=`(beta = Z / sqrt(N), se = 1 / sqrt(N))]

  # Attach chromosome/position and orient alleles in alphabetical order.
  dat <- merge(dat, read_bim(prefix), by = "SNP")
  flip <- dat$A1 > dat$A2
  dat[flip, `:=`(
    beta = -beta,
    Z = -Z,
    A1 = A2,
    A2 = A1
  )]

  dat[, rsid := paste0(chr, ":", position, "_", A1, "_", A2)]
  dat[, .(
    SNP, rsid, chr, position,
    ea = A1, nea = A2, beta, se, p = P, Z, N,
    eaf = NA_real_, pop = prefix,
    id = paste(trait, prefix), trait = trait
  )]
}

# Read each trait pair once and reuse it for both thresholds.
prepare_pair <- function(auxiliary, target, exposure, outcome) {
  pops <- c(auxiliary, target)
  list(
    pops = pops,
    exposure = setNames(lapply(pops, read_gwas, trait = exposure), pops),
    outcome = setNames(lapply(pops, read_gwas, trait = outcome), pops)
  )
}

## =========================================================
## 3. Run one CAMERA analysis
## =========================================================

run_camera <- function(dat, auxiliary, target, exposure, outcome, threshold) {
  pops <- dat$pops

  # 1. Select significant exposure SNPs in the auxiliary and target populations;
  # clump each set using its population-specific LD reference.
  hits_by_pop <- lapply(pops, function(pop) {
    sig <- dat$exposure[[pop]][p < threshold]
    if (nrow(sig) == 0) return(NULL)

    if (nrow(sig) == 1) return(sig)
    clumped <- ieugwasr::ld_clump(
      data.frame(rsid = sig$SNP, pval = sig$p),
      clump_kb = clump_kb,
      clump_r2 = clump_r2,
      clump_p = threshold,
      bfile = bfiles[[ancestry_of(pop)]],
      plink_bin = plink_bin
    )
    sig[SNP %in% clumped$rsid]
  })
  hits_by_pop <- Filter(Negate(is.null), hits_by_pop)
  if (length(hits_by_pop) == 0) stop("No significant exposure SNP")

  # 2. Pool the population-specific clumped top hits for CAMERA;
  # the pooled set defines the regional analysis candidates.
  pooled_hits <- unique(rbindlist(hits_by_pop)[order(p)], by = "SNP")
  if (nrow(pooled_hits) == 0) stop("No IV remains after population-specific clumping")

  # 3. Define +/-250 kb regions around pooled IVs and merge overlapping regions.
  loci <- reduce(GRanges(
    seqnames = pooled_hits$chr,
    ranges = IRanges(
      start = pmax(1L, pooled_hits$position - region_radius),
      end = pooled_hits$position + region_radius
    )
  ))

  exposure_regions <- list()
  outcome_regions <- list()
  raw_exposure <- list()
  raw_outcome <- list()

  for (i in seq_along(loci)) {
    chr_i <- as.character(seqnames(loci)[i])
    start_i <- start(loci)[i]
    end_i <- end(loci)[i]

    exp_region <- setNames(lapply(pops, function(pop) {
      dat$exposure[[pop]][as.character(chr) == chr_i & position >= start_i & position <= end_i]
    }), paste(exposure, pops))

    # Regional CAMERA FEMA uses SNPs shared by both exposure GWASs.
    common_snp <- Reduce(intersect, lapply(exp_region, function(x) x$rsid))
    if (length(common_snp) == 0) next
    exp_region <- lapply(exp_region, function(x) x[rsid %in% common_snp][order(rsid)])

    # Extract outcome associations during harmonization, after IV selection.
    out_region <- setNames(lapply(pops, function(pop) {
      dat$outcome[[pop]][
        as.character(chr) == chr_i & position >= start_i & position <= end_i &
          rsid %in% common_snp
      ]
    }), paste(outcome, pops))

    exp_all <- rbindlist(exp_region)
    lead_snp <- exp_all$rsid[which.max(abs(exp_all$beta / exp_all$se))]
    region_name <- paste0("region_", length(exposure_regions) + 1L)

    exposure_regions[[region_name]] <- exp_region
    outcome_regions[[region_name]] <- out_region
    raw_exposure[[region_name]] <- exp_all[rsid == lead_snp]
    raw_outcome[[region_name]] <- rbindlist(out_region)[rsid == lead_snp]
  }
  if (length(exposure_regions) == 0) stop("No common exposure region for CAMeRa")

  # 4. Use fixed-effect meta-analysis within each region
  # to select a regional instrument shared by both populations.
  fit <- CAMERA$new()
  fit$import_from_local(
    instrument_raw = rbindlist(raw_exposure),
    instrument_outcome = rbindlist(raw_outcome, fill = TRUE),
    instrument_regions = exposure_regions,
    instrument_outcome_regions = outcome_regions,
    exposure_ids = paste(exposure, pops),
    outcome_ids = paste(outcome, pops),
    pops = pops
  )
  fit$fema_regional_instruments(method = "fema")

  # 5. Extract regional-IV outcome associations, harmonize and perform MR.
  fit$make_outcome_local(exp = fit$instrument_fema)
  fit$harmonise(exp = fit$instrument_fema)
  mr <- as.data.frame(fit$cross_estimate())

  pop_col <- intersect(c("pops", "population"), names(mr))[1]
  if (is.na(pop_col)) stop("CAMeRa output has no population column")
  target_result <- mr[as.character(mr[[pop_col]]) == target, , drop = FALSE]
  if (nrow(target_result) != 1) stop("No unique target-population estimate")

  # 6. Format the columns returned by CAMERA.
  get_value <- function(candidates) {
    column <- candidates[candidates %in% names(target_result)][1]
    if (is.na(column)) NA_real_ else as.numeric(target_result[[column]][1])
  }

  beta <- get_value(c("Estimate", "beta"))
  se <- get_value(c("Std. Error", "se"))
  p_value <- get_value(c("Pr(>|t|)", "p", "pval"))
  if (!is.finite(beta) || !is.finite(se) || se <= 0) {
    stop("CAMeRa returned an invalid target-population estimate")
  }
  z <- beta / se
  if (!is.finite(p_value)) p_value <- 2 * pnorm(abs(z), lower.tail = FALSE)

  data.table(
    analysis_id = paste(auxiliary, target, exposure, outcome, threshold, sep = "|"),
    method = "CAMeRa",
    auxiliary = auxiliary,
    target = target,
    exposure = exposure,
    outcome = outcome,
    threshold = threshold,
    beta = beta,
    se = se,
    z = z,
    p = p_value,
    CI_low = beta - 1.96 * se,
    CI_high = beta + 1.96 * se,
    nIV = uniqueN(fit$instrument_fema$rsid),
    n_pooled_IV = nrow(pooled_hits),
    n_regions = length(exposure_regions),
    clump_kb = clump_kb,
    clump_r2 = clump_r2,
    status = "ok",
    error = NA_character_
  )
}

## =========================================================
## 4. Run analyses and append results
## =========================================================

make_failure <- function(auxiliary, target, exposure, outcome, threshold, error) {
  data.table(
    analysis_id = paste(auxiliary, target, exposure, outcome, threshold, sep = "|"),
    method = "CAMeRa", auxiliary, target, exposure, outcome, threshold,
    beta = NA_real_, se = NA_real_, z = NA_real_, p = NA_real_,
    CI_low = NA_real_, CI_high = NA_real_, nIV = NA_integer_,
    n_pooled_IV = NA_integer_, n_regions = NA_integer_,
    clump_kb, clump_r2, status = "failed", error = as.character(error)[1]
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
  existing_results <- fread(result_file)

  # Keep completed analyses, but remove previous failures so they are retried
  # after code or input corrections without creating duplicate result rows.
  successful_results <- existing_results[status == "ok"]
  if (nrow(successful_results) < nrow(existing_results)) {
    fwrite(
      successful_results, result_file,
      sep = "\t", quote = FALSE, na = "NA"
    )
  }
  completed <- unique(successful_results$analysis_id)
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
            run_camera(prepared, auxiliary, target, exposure, outcome, threshold),
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

message("CAMeRa finished: ", result_file)
