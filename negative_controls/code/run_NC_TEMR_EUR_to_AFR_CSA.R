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
# Negative-control TEMR: EUR auxiliary population -> AFR/CSA targets
#
# Exposure: all phenocodes in info.csv except 1747_5 and 1727
# Outcome: 1747_5 and 1727
# Target populations: AFR and CSA
# Thresholds: 5e-5 and 5e-8
# Clumping: 10,000 kb, r2 = 0.001, target-population reference panel
#
# Inputs have undergone MAF, INFO, duplicate-SNP and
# palindromic-SNP QC. Check fields, align alleles across GWASs,
# and perform the numerical processing required by the model.
############################################################

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(XMR)
  library(TEMR)
  library(MRAPSS)
})

set.seed(1)

# 1. Set negative-control paths, population combinations, traits and parameters.
nc_dir <- file.path(.xmr_root, "negative_controls")
formatted_dir <- file.path(nc_dir, "formatted_data")
info_file <- file.path(nc_dir, "info.csv")
result_dir <- file.path(.xmr_root, "negative_controls/results")

plink_bin <- file.path(nc_dir, "plink_Linux")
bfiles <- c(
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

result_file <- file.path(result_dir, "NC_TEMR_EUR_to_AFR_CSA_clump10000kb.tsv")

# Define exposure and outcome traits from info.csv.
if (!file.exists(info_file)) stop("Negative-control info file not found: ", info_file)
info <- fread(info_file)
if (!"phenocode" %in% names(info)) stop("info.csv must contain phenocode")
outcomes <- c("1747_5", "1727")
exposures <- setdiff(unique(as.character(info$phenocode)), outcomes)
exposures <- exposures[!is.na(exposures) & nzchar(exposures)]
if (length(exposures) == 0L) stop("No negative-control exposures found")

# 2. Check input directories and target-population clumping panels.
if (!dir.exists(formatted_dir)) stop("Formatted-data directory not found: ", formatted_dir)
if (!file.exists(plink_bin)) stop("PLINK binary not found: ", plink_bin)

required_bfiles <- unlist(lapply(bfiles, paste0, c(".bed", ".bim", ".fam")))
missing_bfiles <- required_bfiles[!file.exists(required_bfiles)]
if (length(missing_bfiles) > 0) {
  stop("Missing target-population reference files: ", paste(missing_bfiles, collapse = ", "))
}

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

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

# 5. Align exposure and outcome GWASs within each population.
prepare_pair <- function(auxiliary, target, exposure, outcome) {
  exposure_target <- read_gwas(target, exposure)
  exposure_auxiliary <- read_gwas(auxiliary, exposure)
  outcome_target <- read_gwas(target, outcome)
  outcome_auxiliary <- read_gwas(auxiliary, outcome)

  exposure_aligned <- XMR::align_data(
    dat1 = as.data.frame(exposure_target),
    dat2 = as.data.frame(exposure_auxiliary),
    trait1.name = paste0(target, "_", exposure),
    trait2.name = paste0(auxiliary, "_", exposure),
    h2.fix.intercept = FALSE,
    LDSC = FALSE
  )$dat

  outcome_aligned <- XMR::align_data(
    dat1 = as.data.frame(outcome_target),
    dat2 = as.data.frame(outcome_auxiliary),
    trait1.name = paste0(target, "_", outcome),
    trait2.name = paste0(auxiliary, "_", outcome),
    h2.fix.intercept = FALSE,
    LDSC = FALSE
  )$dat

  if (ncol(exposure_aligned) != 9 || ncol(outcome_aligned) != 9) {
    stop("Unexpected XMR::align_data output")
  }

  exposure_aligned <- as.data.table(exposure_aligned)
  outcome_aligned <- as.data.table(outcome_aligned)

  setnames(exposure_aligned, c(
    "SNP", "A1", "A2", "b_exp_target", "b_exp_auxiliary",
    "se_exp_target", "se_exp_auxiliary", "p_exp_target", "p_exp_auxiliary"
  ))
  setnames(outcome_aligned, c(
    "SNP", "A1", "A2", "b_out_target", "b_out_auxiliary",
    "se_out_target", "se_out_auxiliary", "p_out_target", "p_out_auxiliary"
  ))

  # 6. Harmonize outcome effect directions to exposure allele A1.
  dat <- merge(
    exposure_aligned, outcome_aligned,
    by = "SNP", suffixes = c("_exp", "_out")
  )

  same <- with(dat,
    (A1_exp == A1_out & A2_exp == A2_out) |
      (A1_exp == complement_allele(A1_out) &
         A2_exp == complement_allele(A2_out))
  )
  reverse <- with(dat,
    (A1_exp == A2_out & A2_exp == A1_out) |
      (A1_exp == complement_allele(A2_out) &
         A2_exp == complement_allele(A1_out))
  )

  keep <- which(!is.na(same) & !is.na(reverse) & (same | reverse))
  if (length(keep) == 0) stop("No SNP remains after exposure/outcome harmonisation")

  flip <- reverse[keep] & !same[keep]
  dat <- dat[keep]
  dat[flip, `:=`(
    b_out_target = -b_out_target,
    b_out_auxiliary = -b_out_auxiliary
  )]

  dat[, `:=`(A1 = A1_exp, A2 = A2_exp)]

  list(
    data = dat[, .(
      SNP, A1, A2,
      b_exp_target, b_exp_auxiliary,
      se_exp_target, se_exp_auxiliary,
      p_exp_target, p_exp_auxiliary,
      b_out_target, b_out_auxiliary,
      se_out_target, se_out_auxiliary,
      p_out_target, p_out_auxiliary
    )],
    target_exposure_N = median(exposure_target$N, na.rm = TRUE),
    auxiliary_exposure_N = median(exposure_auxiliary$N, na.rm = TRUE),
    target_outcome_N = median(outcome_target$N, na.rm = TRUE),
    auxiliary_outcome_N = median(outcome_auxiliary$N, na.rm = TRUE)
  )
}

# 7. Extract the target1 causal effect from TEMR output.
extract_temr_result <- function(fit) {
  target_result <- fit$resultall$target1
  if (is.null(target_result)) stop("TEMR output has no target1 result")

  beta_raw <- target_result$beta
  beta <- if (!is.null(names(beta_raw)) && "bbe1" %in% names(beta_raw)) {
    scalar_or_na(beta_raw[["bbe1"]])
  } else {
    scalar_or_na(beta_raw)
  }

  list(
    beta = beta,
    stat = scalar_or_na(target_result$stat),
    pvalue = scalar_or_na(target_result$pvalue)
  )
}

# 8. Format each successful or failed analysis as one result row.
make_result <- function(auxiliary, target, exposure, outcome, threshold, status,
                        error = NA_character_, estimate = NULL, prepared = NULL,
                        n_candidates = NA_integer_, n_iv = NA_integer_) {
  data.table(
    analysis_id = paste(
      auxiliary, target, exposure, outcome, threshold_label(threshold), sep = "|"
    ),
    method = "TEMR",
    auxiliary_dataset = auxiliary,
    target_dataset = target,
    exposure = exposure,
    outcome = outcome,
    threshold = threshold,
    beta = if (is.null(estimate)) NA_real_ else estimate$beta,
    stat = if (is.null(estimate)) NA_real_ else estimate$stat,
    beta_pvalue = if (is.null(estimate)) NA_real_ else estimate$pvalue,
    n_candidates = n_candidates,
    n_IV = n_iv,
    target_exposure_N = if (is.null(prepared)) NA_real_ else prepared$target_exposure_N,
    auxiliary_exposure_N = if (is.null(prepared)) NA_real_ else prepared$auxiliary_exposure_N,
    target_outcome_N = if (is.null(prepared)) NA_real_ else prepared$target_outcome_N,
    auxiliary_outcome_N = if (is.null(prepared)) NA_real_ else prepared$auxiliary_outcome_N,
    matrix_order = "target1,auxiliary",
    candidate_rule = "target P<threshold OR auxiliary P<threshold",
    clump_priority = "sum of unsigned target and auxiliary exposure Z",
    rho = "NULL",
    optimization_method = "CG",
    clump_kb = clump_kb,
    clump_r2 = clump_r2,
    status = status,
    error = as.character(error)[1]
  )
}

# 9. Select joint candidate IVs, clump and fit TEMR at one threshold.
run_temr <- function(prepared, auxiliary, target, exposure, outcome, threshold) {
  dat <- prepared$data

  candidate <-
    (is.finite(dat$p_exp_target) & dat$p_exp_target < threshold) |
    (is.finite(dat$p_exp_auxiliary) & dat$p_exp_auxiliary < threshold)
  candidates <- dat[which(candidate)]
  if (nrow(candidates) == 0) stop("No target/auxiliary exposure IV candidate")

  # Rank clumping candidates by the sum of unsigned Z scores from both exposure P values.
  p_target <- pmin(pmax(candidates$p_exp_target, .Machine$double.xmin), 1)
  p_auxiliary <- pmin(pmax(candidates$p_exp_auxiliary, .Machine$double.xmin), 1)
  combined_z <-
    qnorm(p_target / 2, lower.tail = FALSE) +
    qnorm(p_auxiliary / 2, lower.tail = FALSE)
  candidates[, p_combined := pmax(
    2 * pnorm(abs(combined_z), lower.tail = FALSE),
    .Machine$double.xmin
  )]

  clumped <- MRAPSS::clump(
    as.data.frame(candidates),
    IV.Threshold = 1,
    SNP_col = "SNP",
    pval_col = "p_combined",
    clump_kb = clump_kb,
    clump_r2 = clump_r2,
    bfile = bfiles[[target]],
    plink_bin = plink_bin
  )

  n_iv <- if (is.null(clumped)) 0L else nrow(clumped)
  if (n_iv < 3L) stop("TEMR requires at least three IVs; found ", n_iv)

  # TEMR expects the target population first and the auxiliary population second.
  clumped <- as.data.frame(clumped)
  beta_x <- as.matrix(clumped[, c("b_exp_target", "b_exp_auxiliary")])
  beta_y <- as.matrix(clumped[, c("b_out_target", "b_out_auxiliary")])
  se_x <- as.matrix(clumped[, c("se_exp_target", "se_exp_auxiliary")])
  se_y <- as.matrix(clumped[, c("se_out_target", "se_out_auxiliary")])

  storage.mode(beta_x) <- "double"
  storage.mode(beta_y) <- "double"
  storage.mode(se_x) <- "double"
  storage.mode(se_y) <- "double"
  colnames(beta_x) <- colnames(beta_y) <- c("target1", "auxiliary")
  colnames(se_x) <- colnames(se_y) <- c("target1", "auxiliary")

  fit <- TEMR::TEMR(
    beta_x, beta_y, se_x, se_y,
    rho = NULL,
    meth = "CG"
  )
  estimate <- extract_temr_result(fit)
  if (!is.finite(estimate$beta)) stop("TEMR returned a non-finite target estimate")

  make_result(
    auxiliary, target, exposure, outcome, threshold,
    status = "ok", estimate = estimate, prepared = prepared,
    n_candidates = nrow(candidates), n_iv = n_iv
  )
}

# 10. Read completed analysis IDs to support resuming the run.
results_started <- file.exists(result_file) && file.info(result_file)$size > 0
completed <- character(0)
if (results_started) {
  old_results <- fread(result_file)
  if (!"analysis_id" %in% names(old_results)) stop("Existing result file has no analysis_id")
  completed <- unique(old_results$analysis_id)
}

# 11. Append each completed threshold analysis to the result table.
append_result <- function(row) {
  fwrite(
    row, result_file, sep = "\t", quote = FALSE, na = "NA",
    append = results_started, col.names = !results_started
  )
  results_started <<- TRUE
}

# 12. Run both population combinations, all trait pairs and both thresholds.
message(
  "TEMR analyses: ",
  nrow(combinations) * length(exposures) * length(outcomes) * length(thresholds)
)
message("Result file: ", result_file)

for (combination_index in seq_len(nrow(combinations))) {
  auxiliary <- combinations$auxiliary[combination_index]
  target <- combinations$target[combination_index]

  for (exposure in exposures) {
    for (outcome in outcomes) {
      message("[", auxiliary, " -> ", target, "] ", exposure, " -> ", outcome)

      preparation_error <- NA_character_
      prepared <- tryCatch(
        prepare_pair(auxiliary, target, exposure, outcome),
        error = function(e) {
          preparation_error <<- conditionMessage(e)
          NULL
        }
      )

      for (threshold in thresholds) {
        analysis_id <- paste(
          auxiliary, target, exposure, outcome, threshold_label(threshold), sep = "|"
        )
        if (analysis_id %in% completed) {
          message("  skip existing ", threshold_label(threshold))
          next
        }

        result <- if (is.null(prepared)) {
          make_result(
            auxiliary, target, exposure, outcome, threshold,
            status = "failed", error = preparation_error
          )
        } else {
          tryCatch(
            run_temr(prepared, auxiliary, target, exposure, outcome, threshold),
            error = function(e) {
              candidate <-
                (is.finite(prepared$data$p_exp_target) &
                   prepared$data$p_exp_target < threshold) |
                (is.finite(prepared$data$p_exp_auxiliary) &
                   prepared$data$p_exp_auxiliary < threshold)

              make_result(
                auxiliary, target, exposure, outcome, threshold,
                status = "failed", error = conditionMessage(e), prepared = prepared,
                n_candidates = sum(candidate)
              )
            }
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
