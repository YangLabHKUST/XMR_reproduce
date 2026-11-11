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

options(stringsAsFactors = FALSE, warn = 1)
args <- commandArgs(TRUE)
allowed <- c("--root", "--stage", "--targets", "--thresholds", "--arms", "--retry-failed")
if (any(startsWith(args, "--") & !args %in% allowed)) stop("Unknown option")
arg <- function(key, default) {
  i <- match(key, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop("Missing value: ", key)
  args[i + 1L]
}
root <- normalizePath(arg("--root", file.path(.xmr_root, "SHBG_T2D")), mustWork = TRUE)
code_dir <- file.path(root, "code")
suppressPackageStartupMessages({library(XMR); library(XMAP); library(MRAPSS); library(data.table)})
data.table::setDTthreads(1L)
source(file.path(code_dir, "estimate_bgparas_stage2_C0.R"))
stage <- arg("--stage", "all")
targets <- strsplit(arg("--targets", "EAS,AFR,CSA"), ",", fixed = TRUE)[[1]]
thresholds <- as.numeric(strsplit(arg("--thresholds", "5e-5,5e-8"), ",", fixed = TRUE)[[1]])
arms <- strsplit(arg("--arms", "stage2_C0"), ",", fixed = TRUE)[[1]]
stopifnot(stage %in% c("preflight", "background", "selection", "fit", "summarize", "all", "smoke"),
          all(targets %in% c("EAS", "AFR", "CSA")), !anyDuplicated(targets),
          all(thresholds %in% c(5e-5, 5e-8)), !anyDuplicated(thresholds),
          all(arms %in% c("original", "fit_C0", "stage2_C0")), !anyDuplicated(arms))
smoke <- stage == "smoke"
if (smoke) { targets <- "AFR"; thresholds <- 5e-8; arms <- c("original", "fit_C0", "stage2_C0") }
work <- file.path(root, "results", if (smoke) "xmr_smoke" else "xmr_pipeline")
dir.create(work, recursive = TRUE, showWarnings = FALSE)
setwd(work)
formatted_dir <- file.path(root, "formatted_data")
reference_dir <- nc_dir <- file.path(.xmr_root, "negative_controls")
source_pop <- "EUR"; exposure_code <- "30830"; outcome_code <- "250.2"
clump_kb <- 10000L; clump_r2 <- .001
plink_bin <- file.path(reference_dir, "plink_Linux")
bfiles <- setNames(file.path(reference_dir, "1kg_pops", paste0(c("EAS", "AFR", "SAS"), "_phase3_hg38_autosomes")), c("EAS", "AFR", "CSA"))

source(file.path(code_dir, "xmr_data_helpers.R"))
source(file.path(code_dir, "background_data_helpers.R"))
hash <- function(x) digest::digest(x, algo = "sha256")
fh <- function(f) digest::digest(file = f, algo = "sha256")
versions <- sapply(c("XMR", "XMAP", "MRAPSS", "data.table", "Matrix", "mvtnorm", "expm"), function(p) as.character(packageVersion(p)))
signature <- hash(list(
  code = sapply(file.path(code_dir, c("run_SHBG_T2D_XMR.R", "estimate_bgparas_stage2_C0.R",
                                     "xmr_data_helpers.R", "background_data_helpers.R")), fh),
  input_files = lapply(c(formatted_dir, reference_dir), function(p) {
    f <- list.files(p, recursive = FALSE, full.names = TRUE)
    info <- file.info(f)
    data.frame(name = basename(f), size = info$size, modified = as.numeric(info$mtime))
  }), versions = versions))
writeLines(capture.output(sessionInfo()), file.path(work, "sessionInfo.txt"))
writeLines(capture.output(estimate_bgparas_stage2_C0), file.path(work, "stage2_C0_function.txt"))

atomic_rds <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp.", Sys.getpid())
  saveRDS(x, tmp)
  if (!file.rename(tmp, path)) stop("Cannot publish ", path)
}
cache_read <- function(path, deps = NULL) {
  if (!file.exists(path)) return(NULL)
  x <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.null(x) || !identical(x$checksum, hash(x$value))) {
    q <- paste0(path, ".corrupt.", format(Sys.time(), "%Y%m%d%H%M%S"))
    if (!file.rename(path, q)) stop("Cannot quarantine ", path)
    return(NULL)
  }
  if (!identical(x$signature, signature) || (!is.null(deps) && !identical(x$deps, deps))) {
    stop("Checkpoint source/input signature mismatch: ", path, "; do not mix runs. Use a new output directory or explicitly archive this run.")
  }
  x$value
}
cache_write <- function(x, path, deps = NULL) {
  atomic_rds(list(signature = signature, deps = deps, checksum = hash(x), value = x), path)
}
need_cache <- function(path) {
  x <- cache_read(path)
  if (is.null(x)) stop("Missing prerequisite: ", path)
  x
}
tab <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp.", Sys.getpid())
  fwrite(as.data.table(x), tmp, sep = "\t", na = "NA")
  if (!file.rename(tmp, path)) stop("Cannot publish ", path)
}
same <- function(a, b, label) {
  if (!isTRUE(all.equal(a, b, tolerance = 1e-12))) stop("Invariant failed: ", label)
}
raw_check <- function(x, pop) {
  traits <- c(paste0("EUR_", exposure_code), paste0(pop, "_", exposure_code), paste0(pop, "_", outcome_code))
  expected <- paste(traits[c(1,1,1,2,2,3)], traits[c(1,2,3,2,3,3)])
  stopifnot(nrow(x) == 6L, !anyDuplicated(paste(x$trait1, x$trait2)),
            setequal(paste(x$trait1, x$trait2), expected),
            all(is.finite(as.matrix(x[, c("omega11", "omega12", "omega22", "c1", "c12", "c2")]))))
  x
}
qc <- function(x, pop) {
  names(x)[match(c("omega11", "omega12", "omega22"), names(x))] <- c("Omega11", "Omega12", "Omega22")
  XMR::qc_bgparas(x, exposure_code, outcome_code, "EUR", pop, truncation_corr = .95)
}
raw_matrices <- function(x, pop) {
  traits <- c(paste0("EUR_", exposure_code), paste0(pop,"_",exposure_code), paste0(pop,"_",outcome_code))
  build <- function(col) {
    m <- matrix(NA_real_,3,3,dimnames=list(traits,traits))
    for (i in seq_len(nrow(x))) m[x$trait1[i],x$trait2[i]] <- m[x$trait2[i],x$trait1[i]] <- x[[col]][i]
    m
  }
  list(C=build("c12"),Omega0=build("omega12"))
}
bg_path <- function(pop) file.path(work, pop, "background.rds")
iv_path <- function(pop, threshold) file.path(work, pop, threshold_label(threshold), "selection.rds")
fit_path <- function(pop, threshold, arm) file.path(work, pop, threshold_label(threshold), arm, "result.rds")

preflight <- function() {
  stopifnot(file.access(plink_bin, 1) == 0)
  if (system2(plink_bin, "--version") != 0) stop("PLINK failed")
  for (pop in c("EUR", targets)) {
    for (trait in c(exposure_code, outcome_code)) {
      f <- file.path(formatted_dir, paste0(pop, "_", trait))
      h <- fread(f, nrows = 10L)
      stopifnot(all(c("SNP", "A1", "A2", "Z", "N", "P") %in% names(h)), all(h$N > 0))
    }
  }
  for (pop in targets) {
    stopifnot(all(file.exists(paste0(bfiles[[pop]], c(".bed", ".bim", ".fam")))))
    ld <- load_ldscore_for_bg("EUR", pop)
    stopifnot(all(c(attr(ld, "snp_col"), "EUR", pop, "EUR_EUR", paste0(pop, "_", pop), paste0("EUR_", pop)) %in% names(ld)))
    if (!"allele1" %in% names(ld)) {
      if (pop != "CSA") stop("Missing LD-score allele1 for ",pop)
      # The CSA score table has no allele column; LDSC uses the supplied GWAS
      # Z-score orientation for these inputs.
      stopifnot(identical(which(c("A","G") != NULL), integer()))
      message("CSA LD scores have no allele1 column; using the supplied GWAS Z-score orientation.")
    }
  }
  message("Preflight passed; XMR loaded from ", find.package("XMR"))

}

background <- function(pop) {
  path <- bg_path(pop)
  old <- cache_read(path)
  if (!is.null(old)) { message("Reuse background ", pop); return(invisible(old)) }
  ld <- load_ldscore_for_bg("EUR", pop)
  a <- list(ldscore0 = ld, exposure_code = exposure_code, outcome_code = outcome_code,
            pop_ref = "EUR", pop_target = pop, gwas_dir = formatted_dir, snp_col = attr(ld, "snp_col"))
  message("Unconstrained cross-population background: ", pop)
  original <- raw_check(do.call(XMR::estimate_bgparas,
    c(a, list(cross_intercept_zero = FALSE))), pop)
  if (smoke) {
    off <- raw_check(do.call(estimate_bgparas_stage2_C0, c(a, list(cross_intercept_zero = FALSE))), pop)
    same(original, off, "disabled C0 estimator reproduces original")
  }
  message("Stage 2 C0 background: ", pop)
  constrained <- raw_check(do.call(estimate_bgparas_stage2_C0, a), pop)
  cross <- original$pop1 != original$pop2
  same(original[!cross, ], constrained[!cross, ], "same-population regressions unchanged")
  same(original[, c("omega11", "omega22", "c1", "c2")], constrained[, c("omega11", "omega22", "c1", "c2")], "marginal regressions unchanged")
  stopifnot(all(constrained$c12[cross] == 0))
  mats <- list(original = qc(original, pop), stage2_C0 = qc(constrained, pop))
  mats$fit_C0 <- mats$original
  mats$fit_C0$C[1, 2:3] <- mats$fit_C0$C[2:3, 1] <- 0
  for (arm in c("fit_C0", "stage2_C0")) {
    same(diag(mats[[arm]]$C), diag(mats$original$C), "C diagonals")
    same(mats[[arm]]$C[2:3, 2:3], mats$original$C[2:3, 2:3], "target C block")
    stopifnot(all(mats[[arm]]$C[1, 2:3] == 0), all(mats[[arm]]$C[2:3, 1] == 0))
  }
  stopifnot(identical(mats$fit_C0$Omega0, mats$original$Omega0))
  result <- list(raw_original = original, raw_stage2_C0 = constrained, matrices = mats)
  cache_write(result, path)
  invisible(result)
}

selection <- function(pop, threshold, prepared) {
  path <- iv_path(pop, threshold)
  old <- cache_read(path)
  if (!is.null(old)) { message("Reuse IVs ", pop, " ", threshold); return(invisible(old)) }
  candidates <- subset(prepared$data, pval.exp.pop1 < threshold)
  if (!nrow(candidates)) stop("No EUR-selected candidates")
  candidates <- add_two_population_fema(candidates)
  stopifnot(all(is.finite(candidates$pval.exp.fema2)))
  td <- file.path(dirname(path), "clump")
  dir.create(td, recursive = TRUE, showWarnings = FALSE)
  previous <- getwd(); setwd(td); on.exit(setwd(previous))
  clumped <- MRAPSS::clump(candidates, IV.Threshold = 1, SNP_col = "SNP", pval_col = "pval.exp.fema2",
                         clump_kb = clump_kb, clump_r2 = clump_r2, clump_p = 1,
                         bfile = bfiles[[pop]], plink_bin = plink_bin)
  if (is.null(clumped) || !nrow(clumped)) stop("No IV after clumping")
  clumped <- as.data.frame(clumped)
  stopifnot(!anyDuplicated(clumped$SNP))
  clumped$Threshold <- threshold
  clumped$selection_mode <- "EUR_SELECTED_FEMA_P_CLUMP"
  clumped$clump_priority_mode <- "FEMA_P"
  clumped$effect_heterogeneity_qc <- FALSE
  result <- list(data = clumped, n_candidates = nrow(candidates), n_pairwise_merged = prepared$n_pairwise_merged)
  cache_write(result, path)
  invisible(result)
}

fit_one <- function(pop, threshold, arm) {
  bg <- need_cache(bg_path(pop)); iv <- need_cache(iv_path(pop, threshold))
  mat <- bg$matrices[[arm]]
  deps <- hash(list(iv, mat))
  path <- fit_path(pop, threshold, arm)
  old <- cache_read(path, deps)
  if (!is.null(old) && !("--retry-failed" %in% args && old$row$status != "ok")) {
    message("Reuse fit ", pop, " ", threshold, " ", arm); return(invisible(old))
  }
  dat <- iv$data; tolerance <- if (nrow(dat) < 100) 1e-8 else 1e-7
  message("Fit ", pop, " ", threshold, " ", arm, "; nIV=", nrow(dat))
  warnings <- character()
  set.seed(20261018L)
  fit <- tryCatch(withCallingHandlers(XMR::fit_XMR(data = dat, C = mat$C, Omega0 = mat$Omega0,
                  Threshold = threshold, tol1 = tolerance, tol2 = tolerance, min_thres = 1e-2),
                  warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }), error = identity)
  row <- empty_result(pop, threshold, "fit_error")
  row[, `:=`(arm = arm, nIV = nrow(dat), n_candidates = iv$n_candidates,
             n_pairwise_merged = iv$n_pairwise_merged, n_clumped_IV = nrow(dat),
             input_checksum = deps, warnings = paste(unique(warnings), collapse = " | "))]
  posterior <- NULL
  if (inherits(fit, "error")) {
    row[, error := conditionMessage(fit)]
    fit <- list(error = conditionMessage(fit))
  } else {
    row[, `:=`(beta = scalar_or_na(fit$beta), se = scalar_or_na(fit$beta.se), p = scalar_or_na(fit$beta.pvalue),
               tau.sq = scalar_or_na(fit$tau.sq), pi0 = scalar_or_na(fit$pi0), nvalid = scalar_or_na(fit$nvalid),
               IVsignal.sum = scalar_or_na(fit$post$IVsignal.sum),
               sigmax11 = matrix_value(fit$SigmaX, 1, 1), sigmax12 = matrix_value(fit$SigmaX, 1, 2),
               sigmax22 = matrix_value(fit$SigmaX, 2, 2))]
    pi <- as.numeric(fit$post$Pi)
    stopifnot(length(pi) == nrow(dat), all(is.finite(pi)), all(pi >= 0 & pi <= 1))
    same(as.numeric(fit$nvalid), as.numeric(fit$pi0 * nrow(dat)), "nvalid definition")
    posterior <- cbind(dat, validity_posterior = pi)
    if (all(is.finite(c(row$beta, row$se, row$p, row$pi0))) && row$se > 0 && row$p >= 0 && row$p <= 1) {
      row[, `:=`(status = "ok", error = NA_character_)]
    } else row[, error := "Non-finite or invalid fitted summary"]
  }
  out <- list(row = row, fit = fit, posterior = posterior, C = mat$C, Omega0 = mat$Omega0, selection_checksum = hash(dat))
  cache_write(out, path, deps)
  invisible(out)
}

summarize <- function() {
  rows <- list()
  for (pop in targets) {
    bg <- need_cache(bg_path(pop))
    tab(bg$raw_stage2_C0, file.path(work, pop, "background_stage2_C0.tsv"))
    for (threshold in thresholds) {
      iv <- need_cache(iv_path(pop, threshold))
      tab(iv$data, file.path(dirname(iv_path(pop, threshold)), "IVs.tsv"))
      for (arm in arms) {
        x <- need_cache(fit_path(pop, threshold, arm))
        stopifnot(identical(x$selection_checksum, hash(iv$data)))
        tab(x$row, file.path(dirname(fit_path(pop, threshold, arm)), "result.tsv"))
        if (!is.null(x$posterior)) tab(x$posterior, file.path(dirname(fit_path(pop, threshold, arm)), "IV_posterior.tsv"))
        rows[[length(rows) + 1L]] <- x$row
      }
    }
  }
  all_results <- rbindlist(rows, fill = TRUE)
  stopifnot(nrow(all_results) == length(targets) * length(thresholds) * length(arms),
            !anyDuplicated(all_results[, .(population, threshold, arm)]))
  tab(all_results, file.path(work, "results.tsv"))
  if ("stage2_C0" %in% arms) {
    tab(all_results[arm == "stage2_C0"], file.path(root, "results", "crossPopIntercept0",
        "SHBG_T2D_XMR_EUR_selected_FEMA_P_target_clump_noEffectHetQC.tsv"))
  }
  print(all_results[, .(population, threshold, arm, status, beta, se, p, pi0, nIV, nvalid)])
}


if (stage %in% c("preflight", "all", "smoke")) preflight()
if (stage %in% c("background", "all", "smoke")) for (pop in targets) background(pop)
if (stage %in% c("selection", "all", "smoke")) {
  ld_scores <- load_ld_scores()
  for (pop in targets) {
    todo <- thresholds[vapply(thresholds, function(t) is.null(cache_read(iv_path(pop,t))), logical(1))]
    if (length(todo)) {
      # Use the same largest threshold even when running just the strict smoke test.
      requested <- thresholds; thresholds <- c(5e-5,5e-8)
      prepared <- prepare_target(pop); thresholds <- requested
      for (t in todo) selection(pop,t,prepared)
    }
  }
}
if (stage %in% c("fit", "all", "smoke")) for(pop in targets) for(t in thresholds) for(a in arms) fit_one(pop,t,a)
if (stage %in% c("summarize", "all", "smoke")) summarize()
message("Completed stage=",stage,"; output=",work)
