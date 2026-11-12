#!/usr/bin/env Rscript
# Summarize completed simulation fits for Figure 3; no models are rerun.
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || length(args) > 3L) {
  stop("Usage: Rscript prepare_figure3.R RESULTS_ROOT OUTPUT_DIR [AUXILIARY_SIZE_REPLICATES.tsv]")
}
results_root <- args[1]
out_dir <- args[2]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
profiles <- c(Baseline = "baseline", Reduced = "reduced_confounding")
required <- c("method", "threshold", "beta_true", "rho", "n_aux", "n_target",
              "replicate", "status", "beta", "se", "p")
read_results <- function(profile) {
  directory <- file.path(results_root, profiles[[profile]])
  files <- sort(list.files(directory, "^simulation_.*_results[.]tsv$",
                           recursive = TRUE, full.names = TRUE))
  if (!length(files)) stop("No simulation result tables in profile: ", profile)
  rbindlist(lapply(files, function(file) {
    d <- fread(file, na.strings = c("NA", ""), showProgress = FALSE)
    if (!all(required %in% names(d))) stop("Missing required columns in ", basename(file))
    d <- d[, ..required]
    d[, profile := profile]
    d
  }))
}
raw <- rbindlist(lapply(names(profiles), read_results))
key <- c("profile", "method", "threshold", "beta_true", "rho", "n_aux", "n_target", "replicate")
if (anyDuplicated(raw, by = key)) stop("Duplicate fits: provide one result set per configuration.")
valid_fit <- function(d) {
  !is.na(d$status) & d$status == "ok" & is.finite(d$beta) &
    is.finite(d$se) & d$se > 0 & is.finite(d$p) & d$p >= 0 & d$p <= 1
}
raw[, valid := valid_fit(raw)]
methods <- c("XMR_FEMA", "CAMeRa_FEMA", "MRBEE_SOURCE", "MRBEE_TARGET",
             "MRBEE_TL", "TEMR", "MRAPSS_TARGET")
betas <- c(0, .1, .2, .3, .4)
rows <- data.table(profile = c("Baseline", "Baseline", "Reduced"),
                   n_target = c(6000L, 8000L, 6000L))
primary <- merge(raw[n_aux == 40000 & rho == .8 & beta_true %in% betas], rows,
                 by = c("profile", "n_target"))
groups <- c("profile", "n_target", "beta_true", "method", "threshold")
validity <- primary[, .(n_total = .N, n_valid = sum(valid), n_failed = sum(!valid),
                       n_significant = sum(valid & p < .05, na.rm = TRUE)), by = groups]
expected <- rows[, CJ(beta_true = betas, method = methods, threshold = c(5e-8, 5e-5)),
                 by = .(profile, n_target)]
if (nrow(fsetdiff(expected, validity[, ..groups])) || any(validity$n_total != 50L)) {
  stop("Figure 3 requires 50 replicates for every listed method, setting and threshold.")
}
fwrite(validity, file.path(out_dir, "fit_counts.tsv"), sep = "\t")

# A/B/E/F/I/J: ranks are computed separately for each method and setting.
qq <- primary[beta_true == 0 & valid & !(method == "MRAPSS_TARGET" & threshold == 5e-8)]
qq <- qq[order(p, replicate), .(replicate, p, n_valid = .N,
    expected_log10 = -log10((seq_len(.N) - .5) / .N),
    observed_log10 = -log10(pmax(p, .Machine$double.xmin))), by = groups]
stopifnot(all(qq$n_valid == 50L))
qq[, observed_display := pmin(observed_log10, 3)]
fwrite(qq, file.path(out_dir, "qq.tsv"), sep = "\t")

# C/G/K: target-only MR-APSS pools three rho settings, keeping one auxiliary size.
mr <- merge(raw[method == "MRAPSS_TARGET" & threshold == 5e-5 & n_aux == 40000 &
                  rho %in% c(.2, .5, .8) & beta_true %in% betas], rows,
            by = c("profile", "n_target"))
clusters <- mr[, .(rejections = sum(valid & p < .05, na.rm = TRUE),
                   fits = sum(valid), n_rho = uniqueN(rho)),
               by = .(profile, n_target, beta_true, replicate)]
stopifnot(nrow(clusters) == 750L, all(clusters$fits == 3L), all(clusters$n_rho == 3L))
pooled <- mr[, .(n_total = .N, n_valid = sum(valid), n_failed = sum(!valid),
                 n_significant = sum(valid & p < .05, na.rm = TRUE)),
             by = .(profile, n_target, beta_true)]

# Keep the original bootstrap stream order, including slots for unplotted target sizes.
# This reproduces the published intervals when only the Figure 3 configurations are supplied.
set.seed(20261029)
bootstrap_grid <- CJ(profile = names(profiles), n_target = c(4000L, 6000L, 8000L), beta_true = betas)
intervals <- list()
for (i in seq_len(nrow(bootstrap_grid))) {
  setting <- bootstrap_grid[i]
  indices <- matrix(sample.int(50L, 50L * 20000L, replace = TRUE), nrow = 50L)
  d <- clusters[setting, on = .(profile, n_target, beta_true), nomatch = 0][order(replicate)]
  if (!nrow(d)) next
  stopifnot(nrow(d) == 50L)
  boot <- colSums(matrix(d$rejections[indices], nrow = 50L)) / 150
  ci <- quantile(boot, c(.025, .975), names = FALSE, type = 7)
  intervals[[length(intervals) + 1L]] <- cbind(setting, lower = ci[1], upper = ci[2])
}
pooled <- merge(pooled, rbindlist(intervals), by = c("profile", "n_target", "beta_true"))
pooled[, `:=`(method = "MRAPSS_TARGET", threshold = 5e-5,
              ci_method = "replicate-cluster bootstrap percentile")]
exact_ci <- function(d) {
  ci <- t(mapply(function(k, n) if (n > 0) binom.test(k, n)$conf.int else c(NA, NA),
                 d$n_significant, d$n_valid))
  d[, `:=`(lower = ci[, 1], upper = ci[, 2], ci_method = "exact binomial")]
  d
}
power <- validity[(method == "XMR_FEMA" & threshold == 5e-5) |
                    (profile == "Reduced" & method %in% c("CAMeRa_FEMA", "MRBEE_TARGET") & threshold == 5e-8)]
stopifnot(nrow(power) == 25L, all(power$n_valid == 50L))
power <- rbindlist(list(exact_ci(copy(power)), pooled), use.names = TRUE)
power[, `:=`(power = n_significant / n_valid,
             series = fifelse(method == "XMR_FEMA", "XMR_5e05", method))]
setorder(power, profile, n_target, series, beta_true)
fwrite(power, file.path(out_dir, "power_bars.tsv"), sep = "\t")

# D/H/L: use completed baseline fits, or the supplied per-replicate export.
# The export has already removed MR-APSS copies repeated across auxiliary sizes.
aux <- if (length(args) == 3L) fread(args[3]) else copy(raw[profile == "Baseline"])
stopifnot(all(required %in% names(aux)))
aux <- aux[n_target == 6000 & threshold == 5e-5 & rho %in% c(.2, .5, .8) & beta_true %in% betas &
             ((method == "XMR_FEMA" & n_aux %in% c(20000, 40000, 60000)) |
              (method == "MRAPSS_TARGET" & (is.na(n_aux) | n_aux == 40000)))]
aux[, valid := valid_fit(aux)]
aux_groups <- c("method", "rho", "n_aux", "beta_true")
stopifnot(!anyDuplicated(aux, by = c(aux_groups, "replicate")))
aux_counts <- aux[, .(n_total = .N, n_valid = sum(valid),
                     n_significant = sum(valid & p < .05, na.rm = TRUE)), by = aux_groups]
stopifnot(nrow(aux_counts) == 60L, all(aux_counts$n_total == 50L), all(aux_counts$n_valid == 50L))
reference <- aux_counts[method == "MRAPSS_TARGET", .(n_total = sum(n_total), n_valid = sum(n_valid),
                         n_significant = sum(n_significant)), by = beta_true]
reference[, `:=`(method = "MRAPSS_TARGET", rho = NA_real_, n_aux = NA_integer_)]
curves <- rbindlist(list(aux_counts[method == "XMR_FEMA"], reference), use.names = TRUE)
# These exact binomial intervals reproduce the D/H/L intervals in the published figure.
curves <- exact_ci(curves)
curves[, power := n_significant / n_valid]
fwrite(curves, file.path(out_dir, "power_curves.tsv"), sep = "\t")
message("Prepared Figure 3 tables in ", out_dir)
