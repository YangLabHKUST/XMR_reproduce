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
# Target-only MR-APSS: BBJ_X -> BBJ_Y and TW_X -> TW_Y (TW is TPMI).
# Exports one 5e-05 master clump per pair for both methods and both thresholds.
# No EUR GWAS, EUR causal prior, FEMA-selected IVs or borrowed backgrounds.
# Usage: Rscript run_MRAPSS_BBJ_TPMI.R --check-only
#        Rscript run_MRAPSS_BBJ_TPMI.R

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(data.table))
runner_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(runner_file)) runner_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
shared_iv_code_path <- file.path(dirname(gsub("~+~", " ", runner_file, fixed = TRUE)), "shared_target_IVs_BBJ_TPMI.R")
source(shared_iv_code_path, local = TRUE)

default_config <- function() list(
  work_root = file.path(.xmr_root, "BBJ_TPMI"),
  data_root = file.path(.xmr_root, "BBJ_TPMI"),
  reference_root = file.path(.xmr_root, "negative_controls"),
  targets = c("BBJ", "TPMI"), thresholds = c(5e-5, 5e-8),
  clump_kb = 10000L, clump_r2 = 0.001, min_ivs = 4L,
  bg_z2 = 30, bg_nblocks = 200L, seed = 1L,
  check_only = FALSE, retry_failed = FALSE, rerun = FALSE,
  exposures = NULL, outcomes = NULL
)

configure <- function(args = character()) {
  cfg <- default_config()
  flags <- c("check-only", "retry-failed", "rerun")
  values <- c("work-root", "data-root", "reference-root", "formatted-dir",
              "exposure-meta", "outcome-meta", "ldscore", "bfile", "plink",
              "result-dir", "targets", "thresholds", "exposures", "outcomes", "min-ivs")
  i <- 1L
  while (i <= length(args)) {
    arg <- args[i]
    if (arg == "--help") {
      cat("MR-APSS target-only BBJ/TPMI runner\n",
          "Flags: --check-only --retry-failed --rerun\n",
          "Options: ", paste(paste0("--", values), collapse = " "), "\n",
          "List values are comma-separated. TPMI reads TW_* files.\n",
          "One master clump at 5e-05 per pair; saved 5e-08 subset, also used by MRBEE.\n",
          "Defaults: BBJ,TPMI; 5e-05,5e-08; 10,000 kb; r2=0.001; min-IVs=4.\n", sep = "")
      return(NULL)
    }
    key <- sub("^--", "", arg)
    if (!startsWith(arg, "--") || !key %in% c(flags, values)) stop("Unknown option: ", arg)
    field <- gsub("-", "_", key)
    if (key %in% flags) cfg[[field]] <- TRUE else {
      if (i == length(args) || startsWith(args[i + 1L], "--")) stop("Missing value: ", arg)
      i <- i + 1L
      cfg[[field]] <- args[i]
    }
    i <- i + 1L
  }
  split_value <- function(x) if (is.null(x)) NULL else trimws(unlist(strsplit(x, ",", fixed = TRUE)))
  cfg$targets <- toupper(split_value(cfg$targets))
  cfg$targets[cfg$targets == "TW"] <- "TPMI"
  if (!length(cfg$targets) || any(!cfg$targets %in% c("BBJ", "TPMI")) || anyDuplicated(cfg$targets)) stop("Invalid --targets")
  cfg$thresholds <- as.numeric(split_value(as.character(cfg$thresholds)))
  check_shared_thresholds(cfg$thresholds)
  cfg$exposures <- split_value(cfg$exposures)
  cfg$outcomes <- split_value(cfg$outcomes)
  cfg$min_ivs <- as.numeric(cfg$min_ivs)
  if (length(cfg$min_ivs) != 1L || !is.finite(cfg$min_ivs) || cfg$min_ivs < 4 || cfg$min_ivs != as.integer(cfg$min_ivs)) stop("--min-ivs must be an integer >= 4")
  cfg$min_ivs <- as.integer(cfg$min_ivs)
  defaults <- list(
    formatted_dir = file.path(cfg$data_root, "formatted_data"),
    exposure_meta = file.path(cfg$data_root, "raw_data/BBJ/EUR_BBJ_exposure_description.csv"),
    outcome_meta = file.path(cfg$data_root, "raw_data/BBJ/EUR_BBJ_outcome_description.csv"),
    ldscore = file.path(cfg$reference_root, "LDscore_eas_brit_afr.txt"),
    bfile = file.path(cfg$reference_root, "1kg_pops/EAS_phase3_hg38_autosomes"),
    plink = file.path(cfg$reference_root, "plink_Linux"),
    result_dir = file.path(cfg$work_root, "results/sharedIV_masterClump5e-05/MRAPSS")
  )
  for (nm in names(defaults)) if (is.null(cfg[[nm]])) cfg[[nm]] <- defaults[[nm]]
  cfg
}

threshold_label <- function(p) format(p, scientific = TRUE, digits = 15, trim = TRUE)
dataset_code <- function(target) if (target == "TPMI") "TW" else target
analysis_id <- function(target, exposure, outcome, threshold) paste(target, exposure, outcome, threshold_label(threshold), sep = "|")
complement <- function(x) chartr("ATCG", "TAGC", x)

resolve_gwas <- function(cfg, target, trait) {
  stem <- file.path(cfg$formatted_dir, paste0(dataset_code(target), "_", trait))
  hits <- paste0(stem, c("", ".tsv", ".txt", ".gz", ".bgz", ".tsv.gz", ".txt.gz"))
  hits <- hits[file.exists(hits)]
  if (length(hits) != 1L) stop("Expected exactly one GWAS for ", stem, "; found ", length(hits))
  normalizePath(hits, mustWork = TRUE)
}

trait_table <- function(path, code_col) {
  dat <- fread(path, colClasses = "character")
  if (!code_col %in% names(dat)) stop("Missing metadata column ", code_col, ": ", path)
  name_col <- if ("UKB.Phenotype" %in% names(dat)) "UKB.Phenotype" else code_col
  ans <- unique(data.table(code = trimws(dat[[code_col]]), label = trimws(dat[[name_col]])))
  ans <- ans[!is.na(code) & nzchar(code)]
  if (anyDuplicated(ans$code)) stop("Conflicting duplicate trait codes in ", path)
  ans
}

build_pairs <- function(cfg) {
  ex <- trait_table(cfg$exposure_meta, "UKB.code")
  ou <- trait_table(cfg$outcome_meta, "phecode_unique")
  if (nrow(ex) != 7L || nrow(ou) != 15L) stop("Expected the existing 7-exposure / 15-outcome study metadata")
  if (!is.null(cfg$exposures)) {
    if (any(!cfg$exposures %in% ex$code)) stop("Unknown --exposures code")
    ex <- ex[code %chin% cfg$exposures]
  }
  if (!is.null(cfg$outcomes)) {
    if (any(!cfg$outcomes %in% ou$code)) stop("Unknown --outcomes code")
    ou <- ou[code %chin% cfg$outcomes]
  }
  pairs <- CJ(target = cfg$targets, exposure = ex$code, outcome = ou$code, sorted = FALSE)
  pairs[, `:=`(exposure_name = ex$label[match(exposure, ex$code)], outcome_name = ou$label[match(outcome, ou$code)])]
  pairs[, exposure_file := mapply(function(t, e) resolve_gwas(cfg, t, e), target, exposure)]
  pairs[, outcome_file := mapply(function(t, o) resolve_gwas(cfg, t, o), target, outcome)]
  pairs
}

preflight <- function(cfg) {
  required <- c(cfg$exposure_meta, cfg$outcome_meta, cfg$ldscore, cfg$plink, paste0(cfg$bfile, c(".bed", ".bim", ".fam")))
  problems <- character()
  if (!dir.exists(cfg$formatted_dir)) problems <- c(problems, paste("Missing GWAS directory:", cfg$formatted_dir))
  if (any(!file.exists(required))) problems <- c(problems, paste("Missing file:", required[!file.exists(required)]))
  if (file.exists(cfg$plink) && file.access(cfg$plink, 1) != 0) problems <- c(problems, "PLINK is not executable")
  for (pkg in c("data.table", "XMAP", "MRAPSS")) if (!requireNamespace(pkg, quietly = TRUE)) problems <- c(problems, paste("Missing R package:", pkg))
  if (length(problems)) stop(paste(c("Preflight failed (no MR fits started):", problems), collapse = "\n"))
  build_pairs(cfg)
}

read_gwas <- function(path) {
  dat <- fread(path)
  cols <- c("SNP", "A1", "A2", "Z", "N", "P")
  if (!all(cols %in% names(dat))) stop("GWAS missing columns: ", path)
  dat <- dat[, ..cols]
  dat[, `:=`(SNP = as.character(SNP), A1 = toupper(A1), A2 = toupper(A2),
             Z = as.numeric(Z), N = as.numeric(N), P = as.numeric(P))]
  if (anyDuplicated(dat$SNP)) stop("Duplicate SNP identifiers; resolve before running: ", path)
  # Model-required finite observations and unambiguous allele orientation only.
  dat[!is.na(SNP) & nzchar(SNP) & A1 %chin% c("A", "C", "G", "T") &
        A2 %chin% c("A", "C", "G", "T") & A1 != A2 & A1 != complement(A2) &
        is.finite(Z) & is.finite(N) & N > 0 & is.finite(P) & P >= 0 & P <= 1]
}

load_ldscore <- function(path) {
  dat <- fread(path)
  cols <- c("rsid", "allele1", "allele2", "EAS")
  if (!all(cols %in% names(dat))) stop("LD-score file requires rsid, allele1, allele2, EAS")
  ans <- dat[, .(SNP = as.character(rsid), LD_A1 = toupper(allele1), LD_A2 = toupper(allele2), L2 = as.numeric(EAS))]
  if (anyDuplicated(ans$SNP)) stop("Duplicate SNPs in LD-score reference")
  ans[!is.na(SNP) & nzchar(SNP) & LD_A1 %chin% c("A", "C", "G", "T") &
        LD_A2 %chin% c("A", "C", "G", "T") & LD_A1 != LD_A2 & LD_A1 != complement(LD_A2) & is.finite(L2)]
}

align_to_ld <- function(dat, ld) {
  ans <- merge(dat, ld, by = "SNP", sort = FALSE)
  same <- (ans$A1 == ans$LD_A1 & ans$A2 == ans$LD_A2) | (complement(ans$A1) == ans$LD_A1 & complement(ans$A2) == ans$LD_A2)
  flip <- (ans$A1 == ans$LD_A2 & ans$A2 == ans$LD_A1) | (complement(ans$A1) == ans$LD_A2 & complement(ans$A2) == ans$LD_A1)
  keep <- which(!is.na(same) & !is.na(flip) & (same | flip))
  ans$Z[which(flip & !same)] <- -ans$Z[which(flip & !same)]
  ans <- ans[keep]
  ans[, `:=`(A1 = LD_A1, A2 = LD_A2)]
  ans[, c("LD_A1", "LD_A2") := NULL]
  ans
}

pair_data <- function(x, y) {
  dat <- merge(x, y, by = "SNP", suffixes = c(".x", ".y"), sort = FALSE)
  if (!nrow(dat)) stop("No aligned target SNPs")
  if (any(dat$A1.x != dat$A1.y | dat$A2.x != dat$A2.y)) stop("Allele alignment mismatch")
  dat[, .(SNP, A1 = A1.x, A2 = A2.x, b.exp = Z.x/sqrt(N.x), b.out = Z.y/sqrt(N.y),
          se.exp = 1/sqrt(N.x), se.out = 1/sqrt(N.y), pval.exp = P.x, pval.out = P.y, L2 = L2.x)]
}

estimate_entry <- function(x, y, cfg, gc_fun) {
  s <- intersect(x$SNP, y$SNP)
  x <- x[match(s, SNP)]; y <- y[match(s, SNP)]
  idx <- which(x$Z^2 < cfg$bg_z2 & y$Z^2 < cfg$bg_z2)
  if (length(s) <= cfg$bg_nblocks || length(idx) <= cfg$bg_nblocks) stop("Too few SNPs for two-stage LDSC")
  w <- 1/pmax(x$L2, 1)
  fit1 <- gc_fun(as.data.frame(x[idx, .(Z, N)]), as.data.frame(y[idx, .(Z, N)]),
                 x$L2[idx], y$L2[idx], x$L2[idx], reg_w1 = w[idx], reg_w2 = w[idx], reg_wx = w[idx],
                 constrain_intercept = FALSE, nblocks = cfg$bg_nblocks)
  ints <- c(fit1$tau1$coefs[1], fit1$tau2$coefs[1], fit1$theta$coefs[1])
  if (any(!is.finite(ints))) stop("Invalid stage-1 LDSC intercepts")
  fit2 <- gc_fun(as.data.frame(x[, .(Z, N)]), as.data.frame(y[, .(Z, N)]),
                 x$L2, y$L2, x$L2, reg_w1 = w, reg_w2 = w, reg_wx = w,
                 constrain_intercept = TRUE, int1 = ints[1], int2 = ints[2], intx = ints[3], nblocks = cfg$bg_nblocks)
  # XMAP regresses Z1*Z2 on sqrt(N1*N2)*LD; slope is already per-SNP Omega.
  ans <- c(C11 = fit2$tau1$coefs[1], C12 = fit2$theta$coefs[1], C22 = fit2$tau2$coefs[1],
           O11 = fit2$tau1$coefs[2], O12 = fit2$theta$coefs[2], O22 = fit2$tau2$coefs[2])
  if (length(ans) != 6L || any(!is.finite(ans))) stop("Invalid stage-2 LDSC estimates")
  ans
}

estimate_background <- function(x, y, pair, cfg, gc_fun, self_cache) {
  diagonal <- function(dat, label) {
    if (!exists(label, self_cache, inherits = FALSE)) assign(label, estimate_entry(dat, dat, cfg, gc_fun), self_cache)
    get(label, self_cache, inherits = FALSE)
  }
  xx <- diagonal(x, paste(pair$target, pair$exposure))
  yy <- diagonal(y, paste(pair$target, pair$outcome))
  xy <- estimate_entry(x, y, cfg, gc_fun)
  C <- matrix(c(xx["C11"], xy["C12"], xy["C12"], yy["C11"]), 2)
  Omega <- matrix(c(xx["O11"], xy["O12"], xy["O12"], yy["O11"]), 2)
  list(C = C, Omega = Omega, xx = xx, xy = xy, yy = yy,
       genetic_correlation = Omega[1, 2]/sqrt(Omega[1, 1]*Omega[2, 2]),
       adjusted = FALSE)
}

validate_background <- function(bg) {
  for (nm in c("C", "Omega")) {
    mat <- bg[[nm]]
    if (!identical(dim(mat), c(2L, 2L)) || any(!is.finite(mat)) || any(diag(mat) <= 0)) stop("Invalid ", nm, " background")
    if (!isSymmetric(mat) || inherits(try(chol(mat), silent = TRUE), "try-error")) stop(nm, " is not positive definite; raw background saved, no silent clipping")
  }
  invisible(TRUE)
}

fit_threshold <- function(selection, bg, cfg, threshold, fit_fun, prep_error = NA_character_) {
  selected <- selection_at_threshold(selection, threshold)
  state <- list(n_candidates = selected$n_candidates, n_IV = selected$n_IV,
                threshold_effective = selected$threshold_effective, fit = NULL,
                status = "failed", error = NA_character_, warnings = selected$warnings)
  tryCatch(withCallingHandlers({
    iv <- selected$ivs
    validate_shared_iv(iv, threshold)
    if (!is.na(prep_error)) stop(prep_error)
    if (nrow(iv) < cfg$min_ivs) stop("Insufficient IVs: ", nrow(iv), " < configured minimum ", cfg$min_ivs)
    validate_background(bg)
    # Preserve clump's adjusted Threshold. Never overwrite it with the nominal P cutoff.
    state$fit <- fit_fun(as.data.frame(iv), C = bg$C, Omega = bg$Omega, Cor.SelectionBias = TRUE)
    beta <- as.numeric(state$fit$beta); se <- as.numeric(state$fit$beta.se); p <- as.numeric(state$fit$pvalue)
    if (length(beta) != 1L || length(se) != 1L || length(p) != 1L || !is.finite(beta) || !is.finite(se) || se <= 0 || !is.finite(p) || p < 0 || p > 1) stop("Invalid MR-APSS beta/SE/P result")
    state$status <- "ok"
  }, warning = function(w) {
    state$warnings <<- c(state$warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  }), error = function(e) state$error <<- conditionMessage(e))
  state
}

atomic_table <- function(dat, path) {
  tmp <- tempfile(".mrapss-", tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  fwrite(dat, tmp, sep = "\t", quote = TRUE, na = "NA")
  if (!file.rename(tmp, path)) stop("Cannot atomically replace result table: ", path)
}

upsert_result <- function(existing, row) {
  if (!is.null(existing) && nrow(existing)) existing <- existing[!analysis_id %chin% row$analysis_id]
  ans <- rbindlist(list(existing, row), use.names = TRUE, fill = TRUE)
  if (anyDuplicated(ans$analysis_id)) stop("Duplicate result analysis_id")
  ans
}

is_completed <- function(existing, id, retry_failed = FALSE, rerun = FALSE) {
  if (rerun || is.null(existing) || !nrow(existing)) return(FALSE)
  hit <- existing[analysis_id == id]
  if (!nrow(hit)) return(FALSE)
  if (nrow(hit) != 1L) stop("Duplicate stored analysis_id")
  ok <- isTRUE(hit$status == "ok" && is.finite(hit$beta) && is.finite(hit$beta_se) && hit$beta_se > 0 && is.finite(hit$beta_pvalue) && hit$beta_pvalue >= 0 && hit$beta_pvalue <= 1)
  ok || (!retry_failed && hit$status != "ok")
}

scalar <- function(x) if (length(x) == 1L) suppressWarnings(as.numeric(x)) else NA_real_

package_fingerprints <- function() {
  pkgs <- c("MRAPSS", "XMAP", "data.table", "ieugwasr")
  setNames(lapply(pkgs, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) return("unavailable (test harness only)")
    root <- find.package(pkg)
    files <- c(file.path(root, "DESCRIPTION"), list.files(file.path(root, "R"), full.names = TRUE),
               list.files(file.path(root, "libs"), recursive = TRUE, full.names = TRUE))
    files <- files[file.exists(files) & !dir.exists(files)]
    hashes <- tools::md5sum(files)
    names(hashes) <- substring(names(hashes), nchar(root) + 2L)
    list(version = as.character(packageVersion(pkg)), code_md5 = hashes)
  }), pkgs)
}

make_row <- function(pair, threshold, state, cfg, bg = NULL, bg_file = NA_character_, iv_file = NA_character_, x = NULL, y = NULL) {
  # A failed fit never carries a usable beta or P; diagnostics remain available in the RDS.
  fit <- if (state$status == "ok") state$fit else NULL
  bv <- function(name, i, j) if (is.null(bg)) NA_real_ else bg[[name]][i, j]
  data.table(analysis_id = analysis_id(pair$target, pair$exposure, pair$outcome, threshold), method = "MRAPSS",
             target = pair$target, target_dataset = dataset_code(pair$target), population = "EAS", auxiliary = "none",
             exposure = pair$exposure, outcome = pair$outcome, exposure_name = pair$exposure_name, outcome_name = pair$outcome_name,
             exposure_file = pair$exposure_file, outcome_file = pair$outcome_file,
             threshold = threshold, threshold_nominal = threshold, threshold_effective = state$threshold_effective,
             beta = scalar(fit$beta), beta_se = scalar(fit$beta.se), beta_pvalue = scalar(fit$pvalue),
             se = scalar(fit$beta.se), p = scalar(fit$pvalue),
             n_candidates = state$n_candidates, n_IV = state$n_IV, nIV = state$n_IV,
             pi0 = scalar(fit$pi0), expected_foreground_IVs = scalar(fit$pi0)*state$n_IV,
             IVsignal_sum = scalar(fit$post$IVsignal.sum), tau_sq = scalar(fit$tau.sq), sigma_sq = scalar(fit$sigma.sq),
             exposure_N = if (is.null(x)) NA_real_ else median(x$N), outcome_N = if (is.null(y)) NA_real_ else median(y$N),
             C11 = bv("C",1,1), C12 = bv("C",1,2), C22 = bv("C",2,2),
             Omega11 = bv("Omega",1,1), Omega12 = bv("Omega",1,2), Omega22 = bv("Omega",2,2),
             Cor_SelectionBias = TRUE, clump_kb = cfg$clump_kb, clump_r2 = cfg$clump_r2,
             iv_selection = shared_iv_strategy, master_threshold = 5e-5,
             shared_iv_file = iv_file, shared_iv_md5 = if (!is.na(iv_file)) unname(tools::md5sum(iv_file)) else NA_character_,
             effect_scale = "Z/sqrt(N); SE=1/sqrt(N)",
             background_file = bg_file, iv_file = iv_file, min_ivs = cfg$min_ivs,
             status = state$status, error = state$error, warnings = paste(unique(state$warnings), collapse = " | "),
             finished_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
}

run_target <- function(cfg, target, pairs, ld, clump_fun, fit_fun, gc_fun) {
  target_dir <- file.path(cfg$result_dir, target)
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  lock <- file.path(target_dir, ".running.lock")
  if (!dir.create(lock, showWarnings = FALSE)) stop("Existing run lock: ", lock, "; check running processes before removing a stale lock")
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  writeLines(paste("pid", Sys.getpid(), "started", Sys.time()), file.path(lock, "owner.txt"))
  target_value <- target
  pairs <- pairs[pairs$target == target_value]
  if (!nrow(pairs) || any(pairs$target != target_value)) stop("Target routing failure")
  result_file <- file.path(target_dir, paste0("MRAPSS_", target, "_results.tsv"))
  manifest_file <- file.path(target_dir, "run_manifest.rds")
  # Hash all 22 target GWASs even for subset runs, so smoke/full runs share provenance.
  all_cfg <- cfg; all_cfg$targets <- target; all_cfg$exposures <- NULL; all_cfg$outcomes <- NULL
  all_pairs <- build_pairs(all_cfg)
  paths <- sort(unique(c(all_pairs$exposure_file, all_pairs$outcome_file, cfg$ldscore,
                         cfg$exposure_meta, cfg$outcome_meta, cfg$plink, paste0(cfg$bfile, c(".bed", ".bim", ".fam")))))
  message("Hashing input/reference files for safe resume: ", target)
  manifest <- list(schema = 2L, iv_selection = shared_iv_strategy, input_md5 = tools::md5sum(paths),
                   settings = cfg[c("clump_kb", "clump_r2", "min_ivs", "bg_z2", "bg_nblocks", "seed")],
                   code_md5 = tools::md5sum(c(cfg$script_path, shared_iv_code_path)),
                   R_version = R.version.string,
                   packages = package_fingerprints(),
                   clump_code = paste(deparse(clump_fun), collapse = "\n"),
                   fit_code = paste(deparse(fit_fun), collapse = "\n"),
                   gc_code = paste(deparse(gc_fun), collapse = "\n"))
  if (anyNA(manifest$input_md5)) stop("Unable to hash an input file")
  if (file.exists(manifest_file)) {
    if (!identical(readRDS(manifest_file), manifest)) stop("Inputs/code/settings changed. Use a NEW --result-dir; existing results were not altered.")
  } else {
    if (file.exists(result_file)) stop("Existing results lack run_manifest.rds; use a new result directory")
    saveRDS(manifest, manifest_file)
    fwrite(data.table(path = names(manifest$input_md5), md5 = unname(manifest$input_md5)), file.path(target_dir, "input_manifest.tsv"), sep = "\t")
    capture.output(sessionInfo(), file = file.path(target_dir, "sessionInfo.txt"))
  }
  existing <- if (file.exists(result_file)) fread(result_file, colClasses = c(exposure = "character", outcome = "character")) else NULL
  if (!is.null(existing) && (!all(c("analysis_id", "status", "beta", "beta_se", "beta_pvalue", "shared_iv_file", "shared_iv_md5") %in% names(existing)) || anyDuplicated(existing$analysis_id))) stop("Invalid existing result table")
  self_cache <- new.env(parent = emptyenv())
  for (i in seq_len(nrow(pairs))) {
    pair <- pairs[i]
    pending <- cfg$thresholds[!vapply(cfg$thresholds, function(th) is_completed(existing, analysis_id(target, pair$exposure, pair$outcome, th), cfg$retry_failed, cfg$rerun), logical(1))]
    if (!length(pending)) {
      ids <- vapply(cfg$thresholds, function(th) analysis_id(target, pair$exposure, pair$outcome, th), character(1))
      saved <- existing[analysis_id %chin% ids & !is.na(shared_iv_file)]
      if (nrow(saved)) {
        checked <- read_shared_selection(shared_iv_directory(cfg$result_dir, pair), pair, unname(tools::md5sum(manifest_file)))
        for (j in seq_len(nrow(saved))) {
          if (!identical(saved$shared_iv_md5[j], selection_at_threshold(checked, saved$threshold[j])$iv_md5)) stop("Completed result has different shared IVs")
        }
      }
      next
    }
    message("[", target, " ", i, "/", nrow(pairs), "] ", pair$exposure, " -> ", pair$outcome)
    attempt_dir <- tempfile(paste0(pair$exposure, "_", pair$outcome, "_"), tmpdir = target_dir)
    if (!dir.create(attempt_dir)) stop("Cannot create attempt directory")
    x <- y <- dat <- bg <- selection <- NULL; prep_error <- NA_character_; prep_warnings <- character()
    bg_file <- file.path(attempt_dir, "background.rds")
    selection_dir <- shared_iv_directory(cfg$result_dir, pair)
    manifest_md5 <- unname(tools::md5sum(manifest_file))
    tryCatch(withCallingHandlers({
      if (file.exists(file.path(selection_dir, "selection.rds"))) selection <- read_shared_selection(selection_dir, pair, manifest_md5)
      x <- align_to_ld(read_gwas(pair$exposure_file), ld)
      y <- align_to_ld(read_gwas(pair$outcome_file), ld)
      dat <- pair_data(x, y)
      # Export both IV sets before any LDSC/background validation or MR fit can fail.
      if (is.null(selection)) selection <- prepare_shared_selection(dat, cfg, pair, selection_dir, manifest_md5, clump_fun)
      bg <- estimate_background(x, y, pair, cfg, gc_fun, self_cache)
      saveRDS(bg, bg_file)
      validate_background(bg)
    }, warning = function(w) { prep_warnings <<- c(prep_warnings, conditionMessage(w)); invokeRestart("muffleWarning") }),
    error = function(e) prep_error <<- conditionMessage(e))
    for (th in pending) {
      iv_file <- if (is.null(selection)) NA_character_ else selection_at_threshold(selection, th)$iv_file
      # Stable per-analysis seed, independent of resume order.
      set.seed(as.integer((cfg$seed + sum(utf8ToInt(analysis_id(target, pair$exposure, pair$outcome, th)) * seq_along(utf8ToInt(analysis_id(target, pair$exposure, pair$outcome, th))))) %% 2147483646))
      state <- if (!is.null(selection)) fit_threshold(selection, bg, cfg, th, fit_fun, prep_error) else
        list(status = "failed", error = prep_error, warnings = character(), n_candidates = if (is.null(dat)) NA_integer_ else sum(dat$pval.exp < th), n_IV = NA_integer_, threshold_effective = NA_real_, fit = NULL)
      state$warnings <- c(prep_warnings, state$warnings)
      row <- make_row(pair, th, state, cfg, bg, if (file.exists(bg_file)) bg_file else NA_character_, iv_file, x, y)
      saveRDS(list(result = row, fit = state$fit), file.path(attempt_dir, paste0("fit_", threshold_label(th), ".rds")))
      existing <- upsert_result(existing, row)
      atomic_table(existing, result_file)
      message("  ", threshold_label(th), ": ", row$status, "; IVs=", row$n_IV, if (row$status != "ok") paste0("; ", row$error) else "")
    }
    rm(x, y, dat, bg, selection); invisible(gc(verbose = FALSE))
  }
  if (!is.null(existing)) print(existing[, .(n = .N), by = .(target, threshold, status)])
  invisible(existing)
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  cfg <- configure(args)
  if (is.null(cfg)) return(invisible(NULL))
  script_arg <- grep("^--file=", commandArgs(), value = TRUE)
  # Rscript encodes spaces in --file as ~+~ on some R builds (e.g. iCloud paths).
  cfg$script_path <- if (length(script_arg)) normalizePath(gsub("~+~", " ", sub("^--file=", "", script_arg[1]), fixed = TRUE), mustWork = TRUE) else NULL
  pairs <- preflight(cfg)
  message("Preflight passed: ", nrow(pairs), " target/pair combinations x ", length(cfg$thresholds), " thresholds = ", nrow(pairs)*length(cfg$thresholds), " fits")
  message("Output: ", cfg$result_dir)
  if (cfg$check_only) return(invisible(pairs))
  ld <- load_ldscore(cfg$ldscore)
  dir.create(cfg$result_dir, recursive = TRUE, showWarnings = FALSE)
  for (target in cfg$targets) run_target(cfg, target, pairs, ld, MRAPSS::clump, MRAPSS::MRAPSS, XMAP::estimate_gc)
}

if (sys.nframe() == 0L) main()
