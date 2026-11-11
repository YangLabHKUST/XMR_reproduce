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
# TEMR analysis: SHBG (30830) -> T2D (250.2)
#
# Auxiliary population: EUR
# Target populations: EAS, AFR, CSA
# IV thresholds: 5e-5 and 5e-8
# LD clumping: 10000 kb, r2 = 0.001
#
# Each EUR-target pair is aligned independently, then uses the target-or-
# auxiliary exposure P criterion, the official combined P ranking, and the
# target population's LD panel for clumping.
############################################################

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(MRAPSS)
  library(TEMR)
  library(XMR)
})

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
output_dir <- get_arg(
  args, "--output_dir",
  file.path(.xmr_root, "SHBG_T2D/results/TEMR_results_clump10000kb")
)

exposure_code <- "30830"
outcome_code <- "250.2"
auxiliary_pop <- "EUR"
targets <- c("EAS", "AFR", "CSA")
thresholds <- c(5e-5, 5e-8)
clump_kb <- 10000L
clump_r2 <- 0.001

bfiles <- c(
  EUR = file.path(reference_dir, "1kg_pops/EUR_phase3_hg38_autosomes"),
  EAS = file.path(reference_dir, "1kg_pops/EAS_phase3_hg38_autosomes"),
  AFR = file.path(reference_dir, "1kg_pops/AFR_phase3_hg38_autosomes"),
  CSA = file.path(reference_dir, "1kg_pops/SAS_phase3_hg38_autosomes")
)
plink_bin <- get_arg(args, "--plink_bin", file.path(reference_dir, "plink_Linux"))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
result_file <- file.path(output_dir, "SHBG_T2D_TEMR_all_results.tsv")
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
  dat <- dat[
    !is.na(SNP) & nzchar(SNP) & !duplicated(SNP) &
      A1 %chin% c("A", "C", "G", "T") &
      A2 %chin% c("A", "C", "G", "T") & A1 != A2 &
      is.finite(Z) & is.finite(N) & N > 0 &
      is.finite(P) & P >= 0 & P <= 1
  ]
  as.data.frame(dat)
}

# Cache each exposure GWAS, but align and intersect SNPs separately for every
# EUR-target analysis pair.
analysis_pops <- c(auxiliary_pop, targets)
shbg <- setNames(lapply(analysis_pops, read_gwas, trait = exposure_code), analysis_pops)

empty_result <- function(target, aux, threshold, status, error = NA_character_) {
  data.table(
    method = "TEMR", target_population = target, auxiliary_population = aux,
    exposure = exposure_code, outcome = outcome_code,
    threshold = threshold, beta = NA_real_, se = NA_real_, p = NA_real_,
    stat = NA_real_, nIV = NA_integer_,
    n_pairwise_shbg = NA_integer_,
    candidate_rule = "target P < threshold OR auxiliary P < threshold",
    clump_priority = "official unsigned combined exposure P",
    clump_kb = clump_kb, clump_r2 = clump_r2,
    clump_panel = paste0(target, " target panel"),
    status = status, error = error
  )
}

run_temr <- function(target, aux, threshold) {
  out_target <- read_gwas(target, outcome_code)
  out_aux <- read_gwas(aux, outcome_code)

  exp_pair <- XMR::align_data(
    dat1 = shbg[[target]], dat2 = shbg[[aux]],
    trait1.name = paste0(target, "_SHBG"),
    trait2.name = paste0(aux, "_SHBG"),
    h2.fix.intercept = FALSE, LDSC = FALSE
  )$dat
  colnames(exp_pair) <- c(
    "SNP", "A1", "A2", "b.exp.pop1", "b.exp.pop2",
    "se.exp.pop1", "se.exp.pop2", "pval.exp.pop1", "pval.exp.pop2"
  )
  n_pairwise_shbg <- nrow(exp_pair)

  out_pair <- XMR::align_data(
    dat1 = out_target, dat2 = out_aux,
    trait1.name = paste0(target, "_T2D"),
    trait2.name = paste0(aux, "_T2D"),
    h2.fix.intercept = FALSE, LDSC = FALSE
  )$dat
  colnames(out_pair) <- c(
    "SNP", "A1", "A2", "b.out.pop1", "b.out.pop2",
    "se.out.pop1", "se.out.pop2", "pval.out.pop1", "pval.out.pop2"
  )

  dat <- merge(exp_pair, out_pair, by = "SNP", suffixes = c(".exp", ".out"))
  same <- dat$A1.exp == dat$A1.out & dat$A2.exp == dat$A2.out
  swap <- dat$A1.exp == dat$A2.out & dat$A2.exp == dat$A1.out
  strand_same <- dat$A1.exp == complement(dat$A1.out) &
    dat$A2.exp == complement(dat$A2.out)
  strand_swap <- dat$A1.exp == complement(dat$A2.out) &
    dat$A2.exp == complement(dat$A1.out)
  keep <- same | swap | strand_same | strand_swap
  dat <- dat[keep, , drop = FALSE]
  flip <- !same[keep] & (swap[keep] | strand_swap[keep])
  dat$b.out.pop1[flip] <- -dat$b.out.pop1[flip]
  dat$b.out.pop2[flip] <- -dat$b.out.pop2[flip]

  dat <- dat[
    dat$pval.exp.pop1 < threshold | dat$pval.exp.pop2 < threshold,
    , drop = FALSE
  ]
  if (!nrow(dat)) stop("No exposure SNP met the target-or-auxiliary threshold")

  # Preserve the combined-P expression used by the supplied official example.
  combined_z <- qnorm(1 - pmax(dat$pval.exp.pop1, .Machine$double.xmin) / 2) +
    qnorm(1 - pmax(dat$pval.exp.pop2, .Machine$double.xmin) / 2)
  dat$P.combined <- pmax(
    2 * pnorm(combined_z, lower.tail = FALSE),
    .Machine$double.xmin
  )
  clumped <- MRAPSS::clump(
    dat,
    IV.Threshold = 1,
    SNP_col = "SNP", pval_col = "P.combined",
    clump_kb = clump_kb, clump_r2 = clump_r2,
    bfile = bfiles[[target]], plink_bin = plink_bin
  )
  if (is.null(clumped) || nrow(clumped) < 3) {
    stop("Fewer than three IVs remained after clumping")
  }

  beta_x <- as.matrix(clumped[, c("b.exp.pop1", "b.exp.pop2")])
  beta_y <- as.matrix(clumped[, c("b.out.pop1", "b.out.pop2")])
  se_x <- as.matrix(clumped[, c("se.exp.pop1", "se.exp.pop2")])
  se_y <- as.matrix(clumped[, c("se.out.pop1", "se.out.pop2")])
  storage.mode(beta_x) <- storage.mode(beta_y) <- "double"
  storage.mode(se_x) <- storage.mode(se_y) <- "double"
  colnames(beta_x) <- colnames(beta_y) <- c("target1", "auxiliary")
  colnames(se_x) <- colnames(se_y) <- c("target1", "auxiliary")

  fit <- TEMR::TEMR(beta_x, beta_y, se_x, se_y, rho = NULL, meth = "CG")
  fwrite(
    as.data.table(clumped),
    file.path(
      output_dir,
      paste0(
        "TEMR_", target, "_aux_", aux, "_",
        threshold_label(threshold), "_IVs.tsv"
      )
    ),
    sep = "\t", quote = FALSE, na = "NA"
  )

  ans <- empty_result(target, aux, threshold, "ok")
  ans[, `:=`(
    beta = as.numeric(fit$resultall$target1$beta[["bbe1"]]),
    p = as.numeric(fit$resultall$target1$pvalue),
    stat = as.numeric(fit$resultall$target1$stat),
    nIV = nrow(clumped),
    n_pairwise_shbg = n_pairwise_shbg
  )]
  ans
}

results <- list()
k <- 1L
for (threshold in thresholds) {
  for (target in targets) {
    aux <- auxiliary_pop
    message(
      "TEMR | target ", target, " | auxiliary ", aux,
      " | threshold ", threshold_label(threshold)
    )
    results[[k]] <- tryCatch(
      run_temr(target, aux, threshold),
      error = function(e) empty_result(
        target, aux, threshold, "failed", conditionMessage(e)
      )
    )
    k <- k + 1L
  }
}

final <- rbindlist(results, fill = TRUE)
fwrite(final, result_file, sep = "\t", quote = FALSE, na = "NA")
message("TEMR finished: ", result_file)
