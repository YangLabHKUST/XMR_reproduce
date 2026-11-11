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
# Estimate background parameters for the SHBG -> T2D comparison methods.
# Cross-population intercepts are estimated in stage 1 and retained in stage 2.
#
# Exposure: 30830 = SHBG
# Outcome : 250.2 = T2D
############################################################

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(XMR)
  library(XMAP)
  library(data.table)
})

if (!"cross_intercept_zero" %in% names(formals(XMR::estimate_bgparas))) {
  stop("Install the XMR version specified in environment/README.md.")
}
message("Estimating comparator backgrounds with cross_intercept_zero = FALSE.")

## =========================================================
## 1. Paths and fixed parameters
## =========================================================

base_dir <- file.path(.xmr_root, "SHBG_T2D")
formatted_dir <- file.path(base_dir, "formatted_data")
nc_dir <- file.path(dirname(base_dir), "negative_controls")

out_dir <- file.path(.xmr_root, "SHBG_T2D/results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

exposure_code <- "30830"
outcome_code <- "250.2"

out_rdata <- file.path(out_dir, "SHBG-T2D-bgparas.RData")

## =========================================================
## 2. Helpers
## =========================================================

as_numeric_checked <- function(x, label) {
  y <- suppressWarnings(as.numeric(x))
  if (anyNA(y) && !all(is.na(x))) {
    stop("Column contains non-numeric values after conversion: ", label)
  }
  y
}

coerce_ldscore_numeric <- function(ldscore, ref_pop, target_pop) {
  needed <- unique(c(
    ref_pop,
    target_pop,
    paste0(ref_pop, "_", target_pop),
    paste0(target_pop, "_", ref_pop),
    paste0(ref_pop, "_", ref_pop),
    paste0(target_pop, "_", target_pop)
  ))
  needed <- intersect(needed, names(ldscore))
  for (nm in needed) {
    ldscore[[nm]] <- as_numeric_checked(ldscore[[nm]], paste0("ldscore$", nm))
  }
  ldscore
}

load_ldscore_for_bg <- function(ref_pop, target_pop) {
  if ("CSA" %in% c(ref_pop, target_pop)) {
    ldscore <- read.csv(file.path(nc_dir, "LDscore_SAS_EUR.csv"), sep = "\t", check.names = FALSE)
    snp_col <- "SNP"

    ldscore$EUR_CSA <- ldscore$CSA_EUR
    ldscore$EUR_EUR <- ldscore$EUR
    ldscore$CSA_CSA <- ldscore$CSA
  } else {
    ldscore <- fread(file.path(nc_dir, "LDscore_eas_brit_afr.txt"))
    snp_col <- "rsid"

    idx_amb <- which(ldscore$allele1 == comple(ldscore$allele2))
    if (length(idx_amb) > 0) ldscore <- ldscore[-idx_amb, ]

    ldscore$EUR_AFR <- ldscore$AFR_EUR
    ldscore$EAS_AFR <- ldscore$AFR_EAS
    ldscore$EAS_EUR <- ldscore$EUR_EAS
    ldscore$EUR_EUR <- ldscore$EUR
    ldscore$EAS_EAS <- ldscore$EAS
    ldscore$AFR_AFR <- ldscore$AFR
  }

  ldscore <- as.data.frame(ldscore)

  ldscore <- coerce_ldscore_numeric(ldscore, ref_pop, target_pop)
  attr(ldscore, "snp_col") <- snp_col
  ldscore
}

estimate_one_bg <- function(target_pop, ref_pop) {
  message("\n============================================================")
  message("Estimating bgparas for ", ref_pop, " -> ", target_pop)
  message("============================================================")

  ldscore0 <- load_ldscore_for_bg(ref_pop, target_pop)
  snp_col <- attr(ldscore0, "snp_col")

  bgparas <- XMR::estimate_bgparas(
    ldscore0 = ldscore0,
    exposure_code = exposure_code,
    outcome_code = outcome_code,
    pop_ref = ref_pop,
    pop_target = target_pop,
    snp_col = snp_col,
    gwas_dir = formatted_dir,
    cross_intercept_zero = FALSE
  )

  if (is.null(bgparas) || length(dim(bgparas)) < 2 || nrow(bgparas) != 6L) {
    stop("Expected all six background regressions for ", ref_pop, " -> ", target_pop)
  }

  colnames(bgparas) <- c(
    "trait1",
    "trait2",
    "pop1",
    "pop2",
    "Omega11",
    "Omega12",
    "Omega22",
    "c1",
    "c12",
    "c2"
  )

  qc <- qc_bgparas(
    bgparas = bgparas,
    exposure_code = exposure_code,
    outcome_code = outcome_code,
    pop_ref = ref_pop,
    pop_target = target_pop,
    truncation_corr = 0.95
  )

  list(
    target_pop = target_pop,
    ref_pop = ref_pop,
    bgparas = bgparas,
    Omega0 = qc$Omega0,
    C = qc$C
  )
}

## =========================================================
## 3. Run four populations
## =========================================================

bg_plan <- data.frame(
  target_pop = c("EUR", "EAS", "AFR", "CSA"),
  ref_pop = c("EUR", "EUR", "EUR", "EUR"),
  stringsAsFactors = FALSE
)

bg_results <- vector("list", nrow(bg_plan))
names(bg_results) <- bg_plan$target_pop

for (i in seq_len(nrow(bg_plan))) {
  bg_results[[bg_plan$target_pop[i]]] <- estimate_one_bg(
    target_pop = bg_plan$target_pop[i],
    ref_pop = bg_plan$ref_pop[i]
  )
}

bgparas_raw <- lapply(bg_results, `[[`, "bgparas")
Omega0 <- lapply(bg_results, `[[`, "Omega0")
C <- lapply(bg_results, `[[`, "C")
bg_ref_pop <- setNames(bg_plan$ref_pop, bg_plan$target_pop)

## =========================================================
## 4. Save outputs
## =========================================================

for (pop in names(bgparas_raw)) {
  write.table(
    bgparas_raw[[pop]],
    file.path(out_dir, paste0("SHBG-T2D-bgparas-raw-", pop, ".tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

save(
  Omega0,
  C,
  bgparas_raw,
  bg_ref_pop,
  exposure_code,
  outcome_code,
  bg_plan,
  file = out_rdata
)

cat("\nDone. Background parameters saved to:\n", out_rdata, "\n", sep = "")

for (pop in names(Omega0)){
  cat("\n===== ", pop, " | ref = ", bg_ref_pop[[pop]], " =====\n", sep = "")
  print(Omega0[[pop]])
  print(C[[pop]])
}
