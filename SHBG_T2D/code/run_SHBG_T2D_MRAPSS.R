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
# MRAPSS analysis: SHBG (30830) -> T2D (250.2)
#
# Populations: EUR, EAS, AFR, CSA
# IV thresholds: 5e-5 and 5e-8
# LD clumping: 10000 kb, r2 = 0.001
#
# Each population is analysed independently using its own SHBG and T2D
# summary statistics, with population-specific P-value selection and
# clumping against the matching 1000 Genomes LD panel.
############################################################

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(MRAPSS)
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

scalar_or_na <- function(x) {
  if (is.null(x) || !length(x)) return(NA_real_)
  value <- suppressWarnings(as.numeric(x[1]))
  if (is.finite(value)) value else NA_real_
}

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
  file.path(.xmr_root, "SHBG_T2D/results/MRAPSS_results_clump10000kb")
)

exposure_code <- "30830"
outcome_code <- "250.2"
pops <- c("EUR", "EAS", "AFR", "CSA")
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
result_file <- file.path(output_dir, "SHBG_T2D_MRAPSS_all_results.tsv")
if (file.exists(result_file)) file.remove(result_file)

read_gwas <- function(pop, trait) {
  filename <- file.path(formatted_dir, paste0(pop, "_", trait))
  if (!file.exists(filename)) stop("Missing GWAS file: ", filename)

  dat <- fread(filename)
  required <- c("SNP", "A1", "A2", "Z", "N", "P")
  missing <- setdiff(required, names(dat))
  if (length(missing)) {
    stop(filename, " is missing columns: ", paste(missing, collapse = ", "))
  }

  dat[, `:=`(
    SNP = as.character(SNP),
    A1 = toupper(as.character(A1)),
    A2 = toupper(as.character(A2)),
    Z = as.numeric(Z),
    N = as.numeric(N),
    P = as.numeric(P)
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

# Read each population independently. Exposure/outcome alleles are aligned
# later within the same population by XMR::align_data().
shbg <- setNames(lapply(pops, read_gwas, trait = exposure_code), pops)

load_ld_scores <- function() {
  main <- fread(file.path(reference_dir, "LDscore_eas_brit_afr.txt"))
  csa <- fread(file.path(reference_dir, "LDscore_SAS_EUR.csv"))
  list(main = main, csa = csa)
}

ld_scores <- load_ld_scores()

get_ldscore <- function(pop) {
  if (pop == "CSA") {
    ans <- ld_scores$csa[, .(SNP = as.character(SNP), L2 = as.numeric(CSA))]
  } else {
    ans <- ld_scores$main[, .(SNP = as.character(rsid), L2 = as.numeric(get(pop)))]
  }
  ans <- ans[nzchar(SNP) & is.finite(L2)]
  as.data.frame(unique(ans, by = "SNP"))
}

if (!file.exists(background_file)) stop("Missing background parameter file: ", background_file,
  "; run SHBG_T2D/code/prepare_background.R first.")
load(background_file)
if (!exists("C") || !exists("Omega0")) {
  stop("Background file must contain C and Omega0")
}
C_list <- C
Omega0_list <- Omega0
rm(C, Omega0)

empty_result <- function(pop, threshold, status, error = NA_character_) {
  data.table(
    method = "MRAPSS", population = pop,
    exposure = exposure_code, outcome = outcome_code,
    threshold = threshold, beta = NA_real_, se = NA_real_, p = NA_real_,
    nIV = NA_integer_, nvalid = NA_real_, pi0 = NA_real_,
    n_exposure_snps = nrow(shbg[[pop]]),
    n_harmonised = NA_integer_, n_candidates = NA_integer_,
    clump_kb = clump_kb, clump_r2 = clump_r2,
    clump_panel = pop, status = status, error = error
  )
}

run_one <- function(pop, threshold) {
  outcome <- read_gwas(pop, outcome_code)
  aligned <- XMR::align_data(
    dat1 = shbg[[pop]], dat2 = outcome,
    trait1.name = paste0(pop, "_SHBG"),
    trait2.name = paste0(pop, "_T2D"),
    h2.fix.intercept = FALSE, LDSC = FALSE
  )$dat
  colnames(aligned) <- c(
    "SNP", "A1", "A2", "b.exp", "b.out", "se.exp", "se.out",
    "pval.exp", "pval.out"
  )

  mr_dat <- merge(aligned, get_ldscore(pop), by = "SNP")
  mr_dat$Threshold <- threshold
  n_candidates <- sum(mr_dat$pval.exp < threshold, na.rm = TRUE)

  clumped <- MRAPSS::clump(
    mr_dat,
    IV.Threshold = threshold,
    SNP_col = "SNP", pval_col = "pval.exp",
    clump_kb = clump_kb, clump_r2 = clump_r2,
    bfile = bfiles[[pop]], plink_bin = plink_bin
  )
  n_iv <- if (is.null(clumped)) 0L else nrow(clumped)
  if (n_iv < 5L) stop("MRAPSS requires at least five IVs; found ", n_iv)

  fit_dat <- clumped[, c(
    "SNP", "A1", "A2", "b.exp", "se.exp", "pval.exp",
    "b.out", "se.out", "pval.out", "L2", "Threshold"
  )]
  fit <- MRAPSS::MRAPSS(
    fit_dat,
    C = C_list[[pop]][2:3, 2:3],
    Omega = Omega0_list[[pop]][2:3, 2:3],
    Cor.SelectionBias = TRUE
  )

  fwrite(
    as.data.table(clumped),
    file.path(
      output_dir,
      paste0("MRAPSS_", pop, "_", threshold_label(threshold), "_IVs.tsv")
    ),
    sep = "\t", quote = FALSE, na = "NA"
  )

  ans <- empty_result(pop, threshold, "ok")
  ans[, `:=`(
    beta = scalar_or_na(fit$beta),
    se = scalar_or_na(fit$beta.se),
    p = scalar_or_na(fit$pvalue),
    nIV = nrow(clumped),
    nvalid = nrow(clumped) * scalar_or_na(fit$pi0),
    pi0 = scalar_or_na(fit$pi0),
    n_harmonised = nrow(mr_dat),
    n_candidates = n_candidates
  )]
  ans
}

started <- FALSE
for (threshold in thresholds) {
  for (pop in pops) {
    message("MRAPSS | ", pop, " | threshold ", threshold_label(threshold))
    result <- tryCatch(
      run_one(pop, threshold),
      error = function(e) empty_result(pop, threshold, "failed", conditionMessage(e))
    )
    fwrite(
      result, result_file, sep = "\t", quote = FALSE, na = "NA",
      append = started, col.names = !started
    )
    started <- TRUE
  }
}

message("MRAPSS finished: ", result_file)
