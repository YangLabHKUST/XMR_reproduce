# Shared IV contract: MR-APSS produces these files; MRBEE only reads them.
shared_iv_thresholds <- c(5e-5, 5e-8)
shared_iv_strategy <- "MRAPSS_pair_masterClump5e-05_v1"

check_shared_thresholds <- function(thresholds) {
  if (!length(thresholds) || anyNA(thresholds) || anyDuplicated(thresholds) ||
      any(!thresholds %in% shared_iv_thresholds)) stop("--thresholds must be 5e-05 and/or 5e-08")
  invisible(TRUE)
}

shared_iv_directory <- function(root, pair) {
  file.path(root, pair$target, "shared_IVs", paste(pair$exposure, pair$outcome, sep = "_"))
}

validate_shared_iv <- function(iv, threshold) {
  cols <- c("SNP", "A1", "A2", "b.exp", "b.out", "se.exp", "se.out",
            "pval.exp", "pval.out", "L2", "Threshold")
  if (!all(cols %in% names(iv))) stop("Shared IV file is missing required columns")
  if (anyNA(iv$SNP) || any(!nzchar(iv$SNP)) || anyDuplicated(iv$SNP)) stop("Invalid shared IV SNP IDs")
  for (nm in setdiff(cols, c("SNP", "A1", "A2"))) {
    if (!is.numeric(iv[[nm]]) || any(!is.finite(iv[[nm]]))) stop("Invalid shared IV column: ", nm)
  }
  if (any(iv$se.exp <= 0 | iv$se.out <= 0 | iv$pval.exp < 0 | iv$pval.exp >= threshold |
          iv$pval.out < 0 | iv$pval.out > 1)) stop("Invalid shared IV SE or P value")
  if (anyNA(iv$A1) || anyNA(iv$A2) ||
      any(!iv$A1 %chin% c("A", "C", "G", "T") | !iv$A2 %chin% c("A", "C", "G", "T") |
          iv$A1 == iv$A2 | iv$A1 == chartr("ATCG", "TAGC", iv$A2))) stop("Invalid shared IV alleles")
  if (nrow(iv)) {
    effective <- unique(iv$Threshold)
    if (length(effective) != 1L || effective <= 0 || effective > threshold * (1 + 1e-10)) stop("Invalid shared IV Threshold")
    if (threshold == 5e-8 && effective != threshold) stop("5e-08 subset must use Threshold=5e-08")
  }
  invisible(TRUE)
}

read_shared_selection <- function(directory, pair, source_manifest_md5) {
  meta_file <- file.path(directory, "selection.rds")
  if (!file.exists(meta_file)) stop("Missing saved MR-APSS IV selection; run MR-APSS first: ", meta_file)
  meta <- readRDS(meta_file)
  key <- paste(pair$target, pair$exposure, pair$outcome, sep = "|")
  if (!identical(meta$strategy, shared_iv_strategy) || !identical(meta$pair_key, key) ||
      !identical(meta$thresholds, shared_iv_thresholds) || !identical(meta$master_threshold, 5e-5) ||
      !identical(meta$clump_kb, 10000L) || !identical(meta$clump_r2, .001) ||
      !identical(meta$source_manifest_md5, unname(source_manifest_md5))) stop("Shared IV provenance mismatch: ", meta_file)
  labels <- vapply(shared_iv_thresholds, threshold_label, character(1))
  files <- setNames(file.path(directory, paste0("IVs_", labels, ".tsv")), labels)
  hashes <- unname(tools::md5sum(files))
  if (anyNA(hashes) || !identical(hashes, meta$iv_md5)) stop("Shared IV files missing or changed: ", directory)
  tables <- setNames(lapply(seq_along(files), function(i) {
    iv <- fread(files[i], colClasses = c(SNP = "character", A1 = "character", A2 = "character"))
    # fread infers empty numeric columns as logical.
    if (!nrow(iv)) for (nm in setdiff(names(iv), c("SNP", "A1", "A2"))) set(iv, j = nm, value = numeric())
    validate_shared_iv(iv, shared_iv_thresholds[i])
    iv
  }), labels)
  loose <- copy(tables[[1]][pval.exp < 5e-8]); loose[, Threshold := 5e-8]
  if (!isTRUE(all.equal(as.data.frame(loose), as.data.frame(tables[[2]]), check.attributes = FALSE))) stop("5e-08 IVs are not the saved master subset")
  if (!identical(as.numeric(meta$summary$threshold), shared_iv_thresholds) ||
      any(meta$summary$n_IV != vapply(tables, nrow, integer(1)))) stop("Shared IV count metadata mismatch")
  meta$tables <- tables
  meta$iv_files <- files
  meta$selection_file <- meta_file
  meta
}

prepare_shared_selection <- function(dat, cfg, pair, directory, source_manifest_md5, clump_fun) {
  # A successful selection is immutable and reused even after a failed MR fit.
  if (file.exists(file.path(directory, "selection.rds"))) {
    return(read_shared_selection(directory, pair, source_manifest_md5))
  }
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  if (length(list.files(directory, all.files = FALSE))) stop("Incomplete shared IV export; inspect before rerunning: ", directory)
  candidates <- copy(dat[pval.exp < 5e-5])
  setorder(candidates, pval.exp, SNP)
  warnings <- character()
  if (nrow(candidates)) {
    clumped <- withCallingHandlers(as.data.table(clump_fun(as.data.frame(candidates), IV.Threshold = 5e-5,
      SNP_col = "SNP", pval_col = "pval.exp", clump_kb = cfg$clump_kb, clump_r2 = cfg$clump_r2,
      pop = "EAS", bfile = cfg$bfile, plink_bin = cfg$plink)),
      warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
    if (!"SNP" %in% names(clumped) || anyNA(clumped$SNP) || anyDuplicated(clumped$SNP) ||
        any(!clumped$SNP %chin% candidates$SNP)) stop("Invalid master clump SNP output")
    effective <- if (nrow(clumped)) unique(clumped$Threshold) else NA_real_
    if (nrow(clumped) && (length(effective) != 1L || !is.finite(effective) || effective <= 0 || effective > 5e-5)) stop("Invalid master clump-adjusted Threshold")
    # Retrieve the original aligned values, not any values modified by a clump adapter.
    master <- candidates[SNP %chin% clumped$SNP]
    master[, Threshold := as.numeric(effective)]
  } else {
    master <- copy(candidates); master[, Threshold := numeric()]
  }
  stringent <- copy(master[pval.exp < 5e-8])
  # The stringent subset has its own selection cutoff, not the master's adjusted cutoff.
  stringent[, Threshold := 5e-8]
  tables <- list(master, stringent)
  labels <- vapply(shared_iv_thresholds, threshold_label, character(1))
  files <- file.path(directory, paste0("IVs_", labels, ".tsv"))
  for (i in seq_along(tables)) {
    validate_shared_iv(tables[[i]], shared_iv_thresholds[i])
    fwrite(tables[[i]], files[i], sep = "\t", quote = TRUE, na = "NA")
  }
  summary <- data.table(threshold = shared_iv_thresholds,
    n_candidates = vapply(shared_iv_thresholds, function(th) sum(dat$pval.exp < th), integer(1)),
    n_IV = vapply(tables, nrow, integer(1)),
    threshold_effective = vapply(tables, function(iv) if (nrow(iv)) unique(iv$Threshold) else NA_real_, numeric(1)))
  meta <- list(strategy = shared_iv_strategy, pair_key = paste(pair$target, pair$exposure, pair$outcome, sep = "|"),
    thresholds = shared_iv_thresholds, master_threshold = 5e-5, clump_kb = cfg$clump_kb, clump_r2 = cfg$clump_r2,
    source_manifest_md5 = unname(source_manifest_md5), summary = summary,
    iv_md5 = unname(tools::md5sum(files)), warnings = warnings)
  tmp <- tempfile(".selection-", tmpdir = directory)
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(meta, tmp)
  if (!file.rename(tmp, file.path(directory, "selection.rds"))) stop("Cannot publish shared IV selection")
  read_shared_selection(directory, pair, source_manifest_md5)
}

selection_at_threshold <- function(selection, threshold) {
  label <- threshold_label(threshold)
  i <- match(threshold, selection$summary$threshold)
  if (is.na(i)) stop("Threshold missing from shared selection")
  list(ivs = copy(selection$tables[[label]]), n_candidates = selection$summary$n_candidates[i],
       n_IV = selection$summary$n_IV[i], threshold_effective = selection$summary$threshold_effective[i],
       iv_file = unname(selection$iv_files[label]), iv_md5 = selection$iv_md5[i],
       warnings = selection$warnings)
}
