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
