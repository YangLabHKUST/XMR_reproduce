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
# Negative-control XMR: EUR auxiliary exposure -> AFR/CSA targets
#
# Exposure: all phenocodes in info.csv except 1747_5 and 1727
# Outcome: 1747_5 and 1727
# IV threshold: 5e-5 and 5e-8
# LD clumping: 10000 kb, r2 = 0.001, target AFR/CSA reference panel
#
# Formatted GWAS files have already passed QC. This script only
# performs allele harmonisation, auxiliary-P IV selection, FEMA-P target-LD
# clumping, background-parameter handling, and XMR fitting without heterogeneity filtering
# or an additional selection-bias threshold adjustment.
############################################################

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(XMAP)
  library(XMR)
  library(MRAPSS)
})

set.seed(1)

## 1. Paths and analysis settings
nc_dir <- file.path(.xmr_root, "negative_controls")
formatted_dir <- file.path(nc_dir, "formatted_data")
info_file <- file.path(nc_dir, "info.csv")
result_dir <- file.path(.xmr_root, "negative_controls/results/crossPopIntercept0")

plink_bin <- file.path(nc_dir, "plink_Linux")
bfiles <- c(
  AFR = file.path(nc_dir, "1kg_pops/AFR_phase3_hg38_autosomes"),
  CSA = file.path(nc_dir, "1kg_pops/SAS_phase3_hg38_autosomes")
)

sources <- data.table(
  auxiliary = c("EUR", "EUR"),
  auxiliary_name = c("EUR", "EUR"),
  target = c("AFR", "CSA")
)
thresholds <- c(5e-5, 5e-8)
clump_kb <- 10000
clump_r2 <- 0.001
result_file <- file.path(
  result_dir,
  "NC_XMR_EUR_to_AFR_CSA_FEMA_P_target_clump_noEffectHetQC_clump10000kb.tsv"
)

if (!dir.exists(formatted_dir)) stop("Formatted-data directory not found: ", formatted_dir)
if (!file.exists(info_file)) stop("Negative-control info file not found: ", info_file)
if (!file.exists(plink_bin)) stop("PLINK binary not found: ", plink_bin)
missing_bfiles <- unlist(lapply(bfiles, paste0, c(".bed", ".bim", ".fam")))
missing_bfiles <- missing_bfiles[!file.exists(missing_bfiles)]
if (length(missing_bfiles) > 0L) {
  stop("Missing target-population reference files: ", paste(missing_bfiles, collapse = ", "))
}
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

# Read the negative-control trait pairs.
info <- fread(info_file)
if (!"phenocode" %in% names(info)) stop("info.csv must contain phenocode")
outcomes <- c("1747_5", "1727")
exposures <- setdiff(unique(as.character(info$phenocode)), outcomes)
exposures <- exposures[!is.na(exposures) & nzchar(exposures)]
if (length(exposures) == 0L) stop("No negative-control exposures found")
exposure_table <- data.table(trait = exposures, trait_name = exposures)
outcome_table <- data.table(trait = outcomes, trait_name = outcomes)

threshold_label <- function(x) {
  if (isTRUE(all.equal(x, 5e-5))) return("5e-05")
  if (isTRUE(all.equal(x, 5e-8))) return("5e-08")
  format(x, scientific = TRUE)
}

complement_allele <- function(x) chartr("ATCG", "TAGC", toupper(x))

scalar_or_na <- function(x) {
  if (is.null(x) || length(x) == 0L) return(NA_real_)
  value <- suppressWarnings(as.numeric(x[1]))
  if (is.finite(value)) value else NA_real_
}

## 2. Read formatted GWAS files
read_gwas <- function(prefix, trait) {
  stem <- file.path(formatted_dir, paste0(prefix, "_", trait))
  candidates <- c(stem, paste0(stem, ".tsv"), paste0(stem, ".txt"), paste0(stem, ".gz"))
  filename <- candidates[file.exists(candidates)][1]
  if (is.na(filename)) stop("GWAS file not found: ", stem)

  dat <- fread(filename)
  required <- c("SNP", "A1", "A2", "Z", "N", "P")
  if (!all(required %in% names(dat))) {
    stop(filename, " is missing: ", paste(setdiff(required, names(dat)), collapse = ", "))
  }
  dat[, `:=`(
    SNP = as.character(SNP),
    A1 = toupper(as.character(A1)),
    A2 = toupper(as.character(A2)),
    Z = as.numeric(Z),
    N = as.numeric(N),
    P = as.numeric(P)
  )]
  dat[, .(SNP, A1, A2, Z, N, P)]
}

## 3. Load cross-population LD scores
load_xmr_ldscore <- function(target) {
  if (target == "AFR") {
    ld <- fread(file.path(nc_dir, "LDscore_eas_brit_afr.txt"))
    snp_column <- "rsid"
    cross_column <- if ("EUR_AFR" %in% names(ld)) "EUR_AFR" else "AFR_EUR"
  } else if (target == "CSA") {
    ld <- fread(file.path(nc_dir, "LDscore_SAS_EUR.csv"))
    snp_column <- if ("SNP" %in% names(ld)) "SNP" else "rsid"
    cross_column <- if ("EUR_CSA" %in% names(ld)) "EUR_CSA" else "CSA_EUR"
  } else {
    stop("Unsupported target population: ", target)
  }

  required <- c(snp_column, "EUR", target, cross_column)
  if (!all(required %in% names(ld))) {
    stop("LD-score file is missing: ", paste(setdiff(required, names(ld)), collapse = ", "))
  }
  ans <- ld[, .(
    SNP = as.character(get(snp_column)),
    L2.pop1 = as.numeric(EUR),
    L2.pop2 = as.numeric(get(target)),
    L12 = as.numeric(get(cross_column))
  )]
  unique(
    ans[is.finite(L2.pop1) & is.finite(L2.pop2) & is.finite(L12)],
    by = "SNP"
  )
}


# Re-estimate all six background entries. GWAS Z scores are aligned to the
# same allele within each pair; r1_jk*r2_jk is invariant to a common allele flip.
# Cross-population intercepts are constrained during slope estimation.
assert_crosspop_zero <- function(C) {
  if (!is.matrix(C) || !identical(dim(C), c(3L, 3L)) ||
      any(!is.finite(C)) || any(C[1, 2:3] != 0) || any(C[2:3, 1] != 0)) {
    stop("Expected a finite 3x3 C with zero cross-population entries")
  }
  invisible(TRUE)
}

estimate_background_zero <- function(gwas, ld, traits, filename) {
  ld <- as.data.frame(ld)
  ld <- ld[!duplicated(ld$SNP), ]
  gwas <- lapply(gwas, function(d) {
    d <- as.data.frame(d)
    d$SNP <- as.character(d$SNP)
    d$A1 <- toupper(d$A1); d$A2 <- toupper(d$A2)
    d$Z <- as.numeric(d$Z); d$N <- as.numeric(d$N)
    ok <- !is.na(d$SNP) & nzchar(d$SNP) & is.finite(d$Z) &
      is.finite(d$N) & d$N > 0 &
      d$A1 %in% c("A", "C", "G", "T") &
      d$A2 %in% c("A", "C", "G", "T") & d$A1 != d$A2 &
      d$A1 != chartr("ACGT", "TGCA", d$A2)
    d <- d[ok, ]; d[!duplicated(d$SNP), ]
  })
  C <- Omega0 <- matrix(NA_real_, 3L, 3L, dimnames = list(traits, traits))
  entries <- list()
  score_cols <- c("L2.pop1", "L2.pop2", "L2.pop2")
  for (i in 1:3) for (j in i:3) {
    snps <- Reduce(intersect, list(ld$SNP, gwas[[i]]$SNP, gwas[[j]]$SNP))
    d1 <- gwas[[i]][match(snps, gwas[[i]]$SNP), ]
    d2 <- gwas[[j]][match(snps, gwas[[j]]$SNP), ]
    scores <- ld[match(snps, ld$SNP), ]
    comp <- function(a) chartr("ACGT", "TGCA", a)
    same <- (d1$A1 == d2$A1 & d1$A2 == d2$A2) |
      (d1$A1 == comp(d2$A1) & d1$A2 == comp(d2$A2))
    flip <- (d1$A1 == d2$A2 & d1$A2 == d2$A1) |
      (d1$A1 == comp(d2$A2) & d1$A2 == comp(d2$A1))
    d2$Z[flip] <- -d2$Z[flip]
    cross_pop <- i == 1L && j > 1L
    l1 <- as.numeric(scores[[score_cols[i]]])
    l2 <- as.numeric(scores[[score_cols[j]]])
    lx <- if (cross_pop) as.numeric(scores$L12) else l1
    ok <- (same | flip) & is.finite(l1) & is.finite(l2) & is.finite(lx)
    d1 <- d1[ok, ]; d2 <- d2[ok, ]
    l1 <- l1[ok]; l2 <- l2[ok]; lx <- lx[ok]
    ss1 <- d1[, c("Z", "N"), drop = FALSE]
    ss2 <- d2[, c("Z", "N"), drop = FALSE]
    idx <- which(d1$Z^2 < 30 & d2$Z^2 < 30)
    if (length(idx) <= 200L) stop("Too few LDSC SNPs: ", traits[i], " / ", traits[j])
    w1 <- 1 / pmax(l1, 1); w2 <- 1 / pmax(l2, 1)
    if (cross_pop) {
      # Estimate marginal intercepts without fitting a free cross-population intercept.
      f1 <- XMAP::estimate_gc(ss1[idx, ], ss1[idx, ], l1[idx], l1[idx], l1[idx],
        reg_w1 = w1[idx], reg_w2 = w1[idx], reg_wx = w1[idx],
        constrain_intercept = FALSE)
      f2 <- XMAP::estimate_gc(ss2[idx, ], ss2[idx, ], l2[idx], l2[idx], l2[idx],
        reg_w1 = w2[idx], reg_w2 = w2[idx], reg_wx = w2[idx],
        constrain_intercept = FALSE)
      ints <- c(f1$tau1$coefs[1], f2$tau1$coefs[1], 0)
    } else {
      f <- XMAP::estimate_gc(ss1[idx, ], ss2[idx, ], l1[idx], l2[idx], lx[idx],
        reg_w1 = w1[idx], reg_w2 = w2[idx], reg_wx = sqrt(w1[idx] * w2[idx]),
        constrain_intercept = FALSE)
      ints <- c(f$tau1$coefs[1], f$tau2$coefs[1], f$theta$coefs[1])
    }
    if (any(!is.finite(ints))) stop("Invalid first-stage LDSC intercept")
    fit <- XMAP::estimate_gc(ss1, ss2, l1, l2, lx,
      reg_w1 = w1, reg_w2 = w2, reg_wx = sqrt(w1 * w2),
      constrain_intercept = TRUE, int1 = ints[1], int2 = ints[2], intx = ints[3])
    C[i, j] <- C[j, i] <- fit$theta$coefs[1]
    Omega0[i, j] <- Omega0[j, i] <- fit$theta$coefs[2]
    entries[[length(entries) + 1L]] <- data.frame(
      trait1 = traits[i], trait2 = traits[j], cross_population = cross_pop,
      omega = fit$theta$coefs[2], intercept = fit$theta$coefs[1],
      n_stage1 = length(idx), n_stage2 = nrow(ss1))
  }
  assert_crosspop_zero(C)
  if (any(!is.finite(Omega0)) || any(diag(Omega0) <= 0)) stop("Invalid XMR Omega0")
  if (any(diag(C) <= 0) || min(eigen(C, symmetric = TRUE, only.values = TRUE)$values) < -1e-8) {
    stop("Invalid XMR C: non-positive diagonal or not positive semidefinite")
  }
  for (i in 1:2) for (j in (i + 1):3) {
    r <- Omega0[i, j] / sqrt(Omega0[i, i] * Omega0[j, j])
    if (abs(r) > 1) Omega0[i, j] <- Omega0[j, i] <-
      0.95 * sign(r) * sqrt(Omega0[i, i] * Omega0[j, j])
  }
  ans <- list(C = C, Omega0 = Omega0, file = filename,
              cross_population_intercept = 0, entries = do.call(rbind, entries))
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  saveRDS(ans, filename)
  data.table::fwrite(ans$entries, paste0(filename, ".tsv"), sep = "\t")
  ans
}

get_background <- function(auxiliary, exposure, outcome, gwas, ld) {
  traits <- c(paste0(auxiliary, "_", exposure), paste0(target, "_", exposure),
              paste0(target, "_", outcome))
  estimate_background_zero(gwas, ld, traits,
    file.path(result_dir, "background_XMR",
              paste(auxiliary, target, exposure, outcome, "background.rds", sep = "_")))
}

## 5. Harmonise EUR exposure with target exposure and outcome
prepare_pair <- function(auxiliary, exposure, outcome) {
  source_exposure <- read_gwas(auxiliary, exposure)
  target_exposure <- read_gwas(target, exposure)
  target_outcome <- read_gwas(target, outcome)
  ldscore_xmr <- load_xmr_ldscore(target)

  target_pair <- XMR::align_data(
    dat1 = as.data.frame(target_exposure),
    dat2 = as.data.frame(target_outcome),
    trait1.name = paste0(target, "_", exposure),
    trait2.name = paste0(target, "_", outcome),
    h2.fix.intercept = FALSE,
    LDSC = FALSE
  )$dat
  if (ncol(target_pair) != 9L) stop("Unexpected target align_data output")
  target_pair <- as.data.table(target_pair)
  setnames(target_pair, c(
    "SNP", "A1", "A2", "b.exp.pop2", "b.out.pop2",
    "se.exp.pop2", "se.out.pop2", "pval.exp.pop2", "pval.out.pop2"
  ))

  exposure_pair <- XMR::align_data(
    dat1 = as.data.frame(target_exposure),
    dat2 = as.data.frame(source_exposure),
    trait1.name = paste0(target, "_", exposure),
    trait2.name = paste0(auxiliary, "_", exposure),
    h2.fix.intercept = FALSE,
    LDSC = FALSE
  )$dat
  if (ncol(exposure_pair) != 9L) stop("Unexpected exposure align_data output")
  exposure_pair <- as.data.table(exposure_pair)
  setnames(exposure_pair, c(
    "SNP", "A1", "A2", "b.exp.pop2", "b.exp.pop1",
    "se.exp.pop2", "se.exp.pop1", "pval.exp.pop2", "pval.exp.pop1"
  ))

  dat <- merge(
    target_pair,
    exposure_pair[, .(SNP, A1, A2, b.exp.pop1, se.exp.pop1, pval.exp.pop1)],
    by = c("SNP", "A1", "A2")
  )
  dat <- merge(dat, ldscore_xmr, by = "SNP")
  reference_candidates <- dat[
    pval.exp.pop1 < max(thresholds)
  ]
  if (nrow(reference_candidates) == 0L) {
    stop("No auxiliary-exposure candidate below the largest IV threshold")
  }

  list(
    data = reference_candidates,
    n_reference_candidates = nrow(reference_candidates),
    background = get_background(auxiliary, exposure, outcome,
      list(source_exposure, target_exposure, target_outcome), ldscore_xmr),
    source_N = median(source_exposure$N, na.rm = TRUE),
    target_N = median(target_exposure$N, na.rm = TRUE)
  )
}

# The auxiliary P value defines eligibility. Signed FEMA P only ranks eligible SNPs for clumping.
add_fema_p <- function(dat) {
  w1 <- 1 / dat$se.exp.pop1^2
  w2 <- 1 / dat$se.exp.pop2^2
  beta <- (w1 * dat$b.exp.pop1 + w2 * dat$b.exp.pop2) / (w1 + w2)
  se <- sqrt(1 / (w1 + w2))
  dat[, pval.exp.fema := pmax(
    2 * pnorm(abs(beta / se), lower.tail = FALSE),
    .Machine$double.xmin
  )]
  dat
}

make_result <- function(auxiliary, auxiliary_name, exposure, exposure_name,
                        outcome, outcome_name, threshold, status,
                        error = NA_character_, fit = NULL, prepared = NULL,
                        n_reference_candidates = NA_integer_,
                        n_candidates = NA_integer_,
                        n_iv = NA_integer_) {
  bg <- if (is.null(prepared)) NULL else prepared$background
  data.table(
    analysis_id = paste(auxiliary, target, exposure, outcome, threshold_label(threshold), sep = "|"),
    method = "XMR",
    auxiliary = auxiliary,
    auxiliary_name = auxiliary_name,
    target = target,
    exposure = exposure,
    exposure_name = exposure_name,
    outcome = outcome,
    outcome_name = outcome_name,
    threshold = threshold,
    beta = if (is.null(fit)) NA_real_ else scalar_or_na(fit$beta),
    se = if (is.null(fit)) NA_real_ else scalar_or_na(fit$beta.se),
    p = if (is.null(fit)) NA_real_ else scalar_or_na(fit$beta.pvalue),
    tau_sq = if (is.null(fit)) NA_real_ else scalar_or_na(fit$tau.sq),
    pi0 = if (is.null(fit)) NA_real_ else scalar_or_na(fit$pi0),
    n_reference_candidates = n_reference_candidates,
    n_candidates = n_candidates,
    n_IV = if (is.null(fit)) n_iv else scalar_or_na(fit$nIV),
    source_exposure_N = if (is.null(prepared)) NA_real_ else prepared$source_N,
    target_exposure_N = if (is.null(prepared)) NA_real_ else prepared$target_N,
    selection_rule = "auxiliary exposure P < threshold",
    dt = 0,
    modified_threshold = threshold,
    clump_priority = "signed FEMA P",
    clump_reference = target,
    clump_kb = clump_kb,
    clump_r2 = clump_r2,
    background_file = if (is.null(bg)) NA_character_ else bg$file,
    status = status,
    error = as.character(error)[1]
  )
}

run_xmr <- function(prepared, auxiliary, auxiliary_name,
                    exposure, exposure_name, outcome, outcome_name, threshold) {
  candidates <- prepared$data[pval.exp.pop1 < threshold]
  n_candidates <- nrow(candidates)
  if (n_candidates == 0L) stop("No auxiliary-exposure IV candidate")

  # FEMA P ranks the auxiliary-P-selected candidates for target-LD clumping.
  candidates <- add_fema_p(candidates)

  clumped <- MRAPSS::clump(
    as.data.frame(candidates),
    IV.Threshold = 1,
    SNP_col = "SNP",
    pval_col = "pval.exp.fema",
    clump_kb = clump_kb,
    clump_r2 = clump_r2,
    clump_p = 1,
    bfile = bfiles[[target]],
    plink_bin = plink_bin
  )
  n_iv <- if (is.null(clumped)) 0L else nrow(clumped)
  if (n_iv == 0L) stop("No IV remains after target-population LD clumping")
  clumped$Threshold <- threshold

  tolerance <- if (n_iv < 100L) 1e-8 else 1e-7
  assert_crosspop_zero(prepared$background$C)
  fit <- XMR::fit_XMR(
    data = clumped,
    C = prepared$background$C,
    Omega0 = prepared$background$Omega0,
    Threshold = threshold,
    tol1 = tolerance,
    tol2 = tolerance,
    min_thres = 1e-2
  )

  make_result(
    auxiliary, auxiliary_name, exposure, exposure_name,
    outcome, outcome_name, threshold, "ok",
    fit = fit, prepared = prepared,
    n_reference_candidates = prepared$n_reference_candidates,
    n_candidates = n_candidates,
    n_iv = n_iv
  )
}

## 6. Run every pair and append each threshold result immediately
results_started <- file.exists(result_file) && file.info(result_file)$size > 0
completed <- character(0)
if (results_started) completed <- unique(fread(result_file)$analysis_id)

append_result <- function(row) {
  fwrite(
    row, result_file, sep = "\t", quote = FALSE, na = "NA",
    append = results_started, col.names = !results_started
  )
  results_started <<- TRUE
}

total <- nrow(sources) * nrow(exposure_table) * nrow(outcome_table) * length(thresholds)
message("XMR analyses: ", total)
message("Result file: ", result_file)

for (source_index in seq_len(nrow(sources))) {
  auxiliary <- sources$auxiliary[source_index]
  auxiliary_name <- sources$auxiliary_name[source_index]
  target <- sources$target[source_index]

  for (exposure_index in seq_len(nrow(exposure_table))) {
    exposure <- exposure_table$trait[exposure_index]
    exposure_name <- exposure_table$trait_name[exposure_index]

    for (outcome_index in seq_len(nrow(outcome_table))) {
      outcome <- outcome_table$trait[outcome_index]
      outcome_name <- outcome_table$trait_name[outcome_index]
      message("[", auxiliary, " -> ", target, "] ", exposure, " -> ", outcome)

      preparation_error <- NA_character_
      prepared <- tryCatch(
        prepare_pair(auxiliary, exposure, outcome),
        error = function(e) {
          preparation_error <<- conditionMessage(e)
          NULL
        }
      )

      for (threshold in thresholds) {
        analysis_id <- paste(
          auxiliary, target, exposure, outcome, threshold_label(threshold), sep = "|"
        )
        if (analysis_id %in% completed) next

        result <- if (is.null(prepared)) {
          make_result(
            auxiliary, auxiliary_name, exposure, exposure_name,
            outcome, outcome_name, threshold, "failed",
            error = preparation_error
          )
        } else {
          tryCatch(
            run_xmr(
              prepared, auxiliary, auxiliary_name,
              exposure, exposure_name, outcome, outcome_name, threshold
            ),
            error = function(e) make_result(
              auxiliary, auxiliary_name, exposure, exposure_name,
              outcome, outcome_name, threshold, "failed",
              error = conditionMessage(e), prepared = prepared,
              n_reference_candidates = prepared$n_reference_candidates,
              n_candidates = sum(prepared$data$pval.exp.pop1 < threshold)
            )
          )
        }

        append_result(result)
        completed <- c(completed, analysis_id)
        message("  ", threshold_label(threshold), " appended [", result$status, "]")
      }

      rm(prepared)
      invisible(gc(verbose = FALSE))
    }
  }
}

message("XMR finished: ", result_file)
