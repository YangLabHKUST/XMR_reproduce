#!/usr/bin/env Rscript

# Step 6 for the current chr2+chr19 simulation:
#   1. estimate the six official XMR background-parameter pairs;
#   2. assemble C and Omega0 for EUR_X, EAS_X, EAS_Y;
#   3. compare C, Omega0, and genetic correlations with simulation truth.

options(stringsAsFactors = FALSE, warn = 1)

DEFAULT_ROOT <- paste0(
  "simulations/data/",
  "baseline"
)
DEFAULT_BETAS <- c(0, 0.2)
DEFAULT_RHOS <- c(0.4, 0.8)
DEFAULT_N_AUX <- c(40000L, 60000L)
DEFAULT_N_TARGET <- 10000L
DEFAULT_N_ANCHORS <- 50L
DEFAULT_REPLICATES <- 50L
Z2_THRESHOLD <- 30
DEFAULT_MASTER_SEED <- 414L
DEFAULT_H2_EUR_Y <- 0.3
DEFAULT_H2_EUR_X <- 0.2
DEFAULT_H2_EAS_X <- 0.2
DEFAULT_H2_EAS_Y <- 0.3
TRAITS <- c("EUR_X", "EAS_X", "EAS_Y")

# C submatrix in the official XMR order EUR_X, EAS_X, EAS_Y.
C_TRUE <- matrix(
  c(
    1.25, 0, 0,
    0, 1.10, 0.15,
    0, 0.15, 1.05
  ),
  nrow = 3L,
  byrow = TRUE,
  dimnames = list(TRAITS, TRAITS)
)
if (!isSymmetric(C_TRUE) ||
    any(!is.finite(C_TRUE)) ||
    min(eigen(C_TRUE, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
  stop("C_TRUE must be a finite, symmetric positive-definite matrix")
}

# Background genetic-correlation submatrix for U1, U2, V2.
make_r_omega_true <- function(rho) {
  matrix(
    c(
      1.0, rho, 0.3,
      rho, 1.0, 0.4,
      0.3, 0.4, 1.0
    ),
    nrow = 3L,
    byrow = TRUE,
    dimnames = list(TRAITS, TRAITS)
  )
}

make_omega_true <- function(
    m_background, rho, h2_eur_x, h2_eas_x, h2_eas_y) {
  h2 <- c(
    EUR_X = h2_eur_x,
    EAS_X = h2_eas_x,
    EAS_Y = h2_eas_y
  )
  scales <- sqrt(h2 / m_background)
  outer(scales, scales) * make_r_omega_true(rho)
}

float_label <- function(value) {
  label <- sprintf("%.12g", value)
  label <- gsub("-", "m", label, fixed = TRUE)
  label <- gsub(".", "p", label, fixed = TRUE)
  gsub("+", "", label, fixed = TRUE)
}

split_cli_values <- function(arguments, option, defaults) {
  values <- character()
  equals_prefix <- paste0(option, "=")
  for (index in seq_along(arguments)) {
    argument <- arguments[[index]]
    if (startsWith(argument, equals_prefix)) {
      values <- c(values, substring(argument, nchar(equals_prefix) + 1L))
    }
    if (identical(argument, option)) {
      cursor <- index + 1L
      while (cursor <= length(arguments) &&
             !startsWith(arguments[[cursor]], "--")) {
        values <- c(values, arguments[[cursor]])
        cursor <- cursor + 1L
      }
    }
  }
  if (length(values) == 0L) return(defaults)
  values <- trimws(unlist(strsplit(values, ",", fixed = TRUE)))
  values[nzchar(values)]
}

parse_numeric_grid <- function(arguments, option, defaults) {
  text <- split_cli_values(arguments, option, as.character(defaults))
  values <- suppressWarnings(as.numeric(text))
  if (length(values) == 0L || any(!is.finite(values))) {
    stop(option, " requires one or more finite numeric values", call. = FALSE)
  }
  unique(values)
}

parse_integer_grid <- function(arguments, option, defaults) {
  values <- parse_numeric_grid(arguments, option, defaults)
  if (any(values <= 0) || any(values != floor(values))) {
    stop(option, " requires positive integer values", call. = FALSE)
  }
  as.integer(values)
}

parse_master_seed <- function(arguments) {
  text <- split_cli_values(arguments, "--seed", as.character(DEFAULT_MASTER_SEED))
  value <- suppressWarnings(as.numeric(text))
  if (length(value) != 1L || !is.finite(value) || value < 0 ||
      value != floor(value) || value > 2147483646) {
    stop("--seed requires one integer between 0 and 2147483646", call. = FALSE)
  }
  as.integer(value)
}

derive_replicate_seed <- function(master_seed, replicate) {
  modulus <- 2147483646
  mixed <- (as.double(master_seed) + 4000012 +
    9176 * as.double(replicate)) %% modulus
  as.integer(mixed + 1)
}

parse_cli_args <- function(arguments = commandArgs(trailingOnly = TRUE)) {
  if (any(arguments %in% c("-h", "--help"))) {
    cat(
      paste0(
        "Usage: Rscript step6_estimate_evaluate_xmr_background.R ",
        "[--root DIR] [--betas 0 0.2] ",
        "[--rhos 0.4 0.8] [--n-aux 40000 60000] ",
        "[--n-target 10000] [--n-anchors 50] ",
        "[--h2-eur-y 0.3] [--h2-eur-x 0.2] ",
        "[--h2-eas-x 0.2] [--h2-eas-y 0.3] ",
        "[--expected-replicates 50] [--seed 414]\n"
      )
    )
    return(NULL)
  }
  betas <- parse_numeric_grid(arguments, "--betas", DEFAULT_BETAS)
  rhos <- parse_numeric_grid(arguments, "--rhos", DEFAULT_RHOS)
  n_aux <- parse_integer_grid(arguments, "--n-aux", DEFAULT_N_AUX)
  n_target <- parse_integer_grid(
    arguments, "--n-target", DEFAULT_N_TARGET
  )
  n_anchors <- parse_integer_grid(
    arguments, "--n-anchors", DEFAULT_N_ANCHORS
  )
  expected_replicates <- parse_integer_grid(
    arguments, "--expected-replicates", DEFAULT_REPLICATES
  )
  master_seed <- parse_master_seed(arguments)
  h2_eur_y <- parse_numeric_grid(
    arguments, "--h2-eur-y", DEFAULT_H2_EUR_Y
  )
  h2_eur_x <- parse_numeric_grid(
    arguments, "--h2-eur-x", DEFAULT_H2_EUR_X
  )
  h2_eas_x <- parse_numeric_grid(
    arguments, "--h2-eas-x", DEFAULT_H2_EAS_X
  )
  h2_eas_y <- parse_numeric_grid(
    arguments, "--h2-eas-y", DEFAULT_H2_EAS_Y
  )
  root <- split_cli_values(
    arguments, "--root", Sys.getenv("SIM_ROOT", unset = DEFAULT_ROOT)
  )
  if (length(root) != 1L || !nzchar(root)) {
    stop("--root requires one directory", call. = FALSE)
  }
  if (length(expected_replicates) != 1L) {
    stop("--expected-replicates requires one integer", call. = FALSE)
  }
  if (length(n_anchors) != 1L) {
    stop("--n-anchors requires one integer", call. = FALSE)
  }
  h2_values <- c(
    h2_eur_y = h2_eur_y,
    h2_eur_x = h2_eur_x,
    h2_eas_x = h2_eas_x,
    h2_eas_y = h2_eas_y
  )
  if (length(h2_values) != 4L || any(h2_values <= 0 | h2_values >= 1)) {
    stop("Each h2 option requires one value in (0, 1)", call. = FALSE)
  }
  if (any(abs(rhos) >= 1)) {
    stop("Every --rhos value must satisfy -1 < rho < 1", call. = FALSE)
  }
  for (rho in rhos) {
    minimum <- min(eigen(make_r_omega_true(rho), symmetric = TRUE)$values)
    if (minimum <= 0) {
      stop(
        sprintf(
          "rho=%g makes the true background correlation matrix non-positive-definite",
          rho
        ),
        call. = FALSE
      )
    }
  }
  list(
    root = normalizePath(root, mustWork = FALSE),
    betas = betas,
    rhos = rhos,
    n_aux = n_aux,
    n_target = n_target,
    n_anchors = n_anchors[[1L]],
    h2_eur_y = h2_eur_y[[1L]],
    h2_eur_x = h2_eur_x[[1L]],
    h2_eas_x = h2_eas_x[[1L]],
    h2_eas_y = h2_eas_y[[1L]],
    master_seed = master_seed,
    expected_replicates = expected_replicates[[1L]],
    overwrite = "--overwrite" %in% arguments
  )
}

setting_id <- function(beta, rho, n_aux, n_target) {
  paste0(
    "beta", float_label(beta), "_rho", float_label(rho),
    "_naux", n_aux, "_ntarget", n_target
  )
}

scenario_paths <- function(
    root, beta, rho, n_aux, n_target, n_anchors,
    h2_eur_y, h2_eur_x, h2_eas_x, h2_eas_y) {
  current_setting <- setting_id(beta, rho, n_aux, n_target)
  h2_label <- paste0(
    "h2eury", float_label(h2_eur_y),
    "_eurx", float_label(h2_eur_x),
    "_easx", float_label(h2_eas_x),
    "_easy", float_label(h2_eas_y)
  )
  list(
    input_glob = file.path(
      root,
      "summary_statistics",
      current_setting,
      sprintf(
        "xmr_sumstats_*_%s_*core%ds*_fixedIV_rep*.tsv.gz",
        h2_label, n_anchors
      )
    ),
    out_dir = file.path(
      root,
      "xmr_background",
      current_setting
    ),
    setting_id = current_setting
  )
}

require_simulation_packages <- function() {
  required <- c("data.table", "ggplot2", "XMAP", "XMR")
  missing <- required[!vapply(
    required, requireNamespace, logical(1), quietly = TRUE
  )]
  if (length(missing) > 0L) {
    stop(
      "Missing R package(s): ", paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
}

parse_replicate <- function(path) {
  fields <- regmatches(
    basename(path),
    regexec("rep([0-9]+)", basename(path), perl = TRUE)
  )[[1L]]
  if (length(fields) != 2L) stop("Cannot parse replicate: ", path)
  as.integer(fields[[2L]])
}

required_columns <- c(
  "SNP", "A1", "A2", "HAS_BACKGROUND",
  "b.exp.pop1", "se.exp.pop1", "N.exp.pop1",
  "b.exp.pop2", "se.exp.pop2", "N.exp.pop2",
  "b.out.pop2", "se.out.pop2", "N.out.pop2",
  "L2.pop1", "L2.pop2", "L12"
)

read_sumstats <- function(path) {
  header <- names(data.table::fread(path, nrows = 0L, showProgress = FALSE))
  missing <- setdiff(required_columns, header)
  if (length(missing) > 0L) {
    stop("Missing columns in ", path, ": ", paste(missing, collapse = ", "))
  }
  dat <- data.table::fread(
    path, select = required_columns, showProgress = FALSE
  )
  dat[, `:=`(
    SNP = as.character(SNP),
    A1 = toupper(as.character(A1)),
    A2 = toupper(as.character(A2)),
    HAS_BACKGROUND = as.integer(HAS_BACKGROUND)
  )]
  if (anyDuplicated(dat$SNP) || anyNA(dat$SNP)) {
    stop("Missing or duplicated SNP IDs in ", path)
  }
  numeric_columns <- setdiff(required_columns, c("SNP", "A1", "A2"))
  for (column in numeric_columns) {
    data.table::set(dat, j = column, value = as.numeric(dat[[column]]))
  }
  if (any(!is.finite(as.matrix(dat[, ..numeric_columns])))) {
    stop("Non-finite values in ", path)
  }
  dat
}

remove_ambiguous <- function(dat) {
  dat[!(paste0(A1, A2) %in% c("AT", "TA", "CG", "GC"))]
}

write_xmr_inputs <- function(dat, directory) {
  trait_columns <- list(
    EUR_X = c("b.exp.pop1", "se.exp.pop1", "N.exp.pop1"),
    EAS_X = c("b.exp.pop2", "se.exp.pop2", "N.exp.pop2"),
    EAS_Y = c("b.out.pop2", "se.out.pop2", "N.out.pop2")
  )
  z_values <- list()
  for (trait in TRAITS) {
    columns <- trait_columns[[trait]]
    z <- dat[[columns[[1L]]]] / dat[[columns[[2L]]]]
    data.table::fwrite(
      data.table::data.table(
        SNP = dat$SNP,
        A1 = dat$A1,
        A2 = dat$A2,
        Z = z,
        N = dat[[columns[[3L]]]]
      ),
      file.path(directory, trait),
      sep = "\t",
      quote = FALSE
    )
    z_values[[trait]] <- z
  }
  z_values
}

make_ldscore <- function(dat) {
  data.frame(
    rsid = dat$SNP,
    allele1 = dat$A1,
    allele2 = dat$A2,
    EUR = dat$L2.pop1,
    EAS = dat$L2.pop2,
    EUR_EUR = dat$L2.pop1,
    EUR_EAS = dat$L12,
    EAS_EUR = dat$L12,
    EAS_EAS = dat$L2.pop2,
    check.names = FALSE
  )
}

# Two-stage LDSC with zero intercepts for EUR_X--EAS_X and EUR_X--EAS_Y.
estimate_bgparas_crosspop_zero <- function(dat, z_values) {
  ns <- list(dat$N.exp.pop1, dat$N.exp.pop2, dat$N.out.pop2)
  ls <- list(dat$L2.pop1, dat$L2.pop2, dat$L2.pop2)
  pops <- c("EUR", "EAS", "EAS")
  rows <- list()
  for (i in 1:3) for (j in i:3) {
    ss1 <- data.frame(Z = z_values[[TRAITS[i]]], N = ns[[i]])
    ss2 <- data.frame(Z = z_values[[TRAITS[j]]], N = ns[[j]])
    idx <- which(ss1$Z^2 < Z2_THRESHOLD & ss2$Z^2 < Z2_THRESHOLD)
    if (length(idx) <= 200L) stop("Too few stage-1 LDSC SNPs: ", TRAITS[i], " / ", TRAITS[j])
    l1 <- ls[[i]]; l2 <- ls[[j]]
    cross <- pops[i] != pops[j]
    lx <- if (cross) dat$L12 else l1
    w1 <- 1 / pmax(l1, 1); w2 <- 1 / pmax(l2, 1)
    if (cross) {
      f1 <- XMAP::estimate_gc(ss1[idx, ], ss1[idx, ], l1[idx], l1[idx], l1[idx],
        reg_w1 = w1[idx], reg_w2 = w1[idx], reg_wx = w1[idx], constrain_intercept = FALSE)
      f2 <- XMAP::estimate_gc(ss2[idx, ], ss2[idx, ], l2[idx], l2[idx], l2[idx],
        reg_w1 = w2[idx], reg_w2 = w2[idx], reg_wx = w2[idx], constrain_intercept = FALSE)
      ints <- c(f1$tau1$coefs[1], f2$tau1$coefs[1], 0)
    } else {
      f <- XMAP::estimate_gc(ss1[idx, ], ss2[idx, ], l1[idx], l2[idx], lx[idx],
        reg_w1 = w1[idx], reg_w2 = w2[idx], reg_wx = sqrt(w1[idx]*w2[idx]),
        constrain_intercept = FALSE)
      ints <- c(f$tau1$coefs[1], f$tau2$coefs[1], f$theta$coefs[1])
    }
    if (any(!is.finite(ints))) stop("Invalid stage-1 LDSC intercept")
    f <- XMAP::estimate_gc(ss1, ss2, l1, l2, lx,
      reg_w1 = w1, reg_w2 = w2, reg_wx = sqrt(w1*w2),
      constrain_intercept = TRUE, int1 = ints[1], int2 = ints[2], intx = ints[3])
    if (cross && f$theta$coefs[1] != 0) stop("Cross-population intercept must be zero")
    rows[[length(rows)+1L]] <- data.frame(
      trait1 = TRAITS[i], trait2 = TRAITS[j], pop1 = pops[i], pop2 = pops[j],
      omega11 = f$tau1$coefs[2], omega12 = f$theta$coefs[2], omega22 = f$tau2$coefs[2],
      c1 = f$tau1$coefs[1], c12 = f$theta$coefs[1], c2 = f$tau2$coefs[1])
  }
  do.call(rbind, rows)
}

standardize_bgparas <- function(bgparas, replicate) {
  expected <- c(
    "trait1", "trait2", "pop1", "pop2",
    "omega11", "omega12", "omega22", "c1", "c12", "c2"
  )
  if (nrow(bgparas) != 6L || ncol(bgparas) != length(expected)) {
    stop("Background estimator did not return six complete pairs")
  }
  result <- data.table::as.data.table(bgparas)
  data.table::setnames(result, expected)
  result[, replicate := replicate]
  data.table::setcolorder(result, c("replicate", expected))
  result
}

add_stage_counts <- function(bgparas, z_values) {
  bgparas[, n_stage1 := mapply(
    function(first, second) {
      z1 <- z_values[[first]]
      z2 <- z_values[[second]]
      sum(z1^2 < Z2_THRESHOLD & z2^2 < Z2_THRESHOLD)
    },
    trait1,
    trait2,
    USE.NAMES = FALSE
  )]
  bgparas[, n_stage2 := length(z_values[[1L]])]
  bgparas
}

matrix_diagnostics <- function(value) {
  eigenvalues <- eigen(
    (value + t(value)) / 2,
    symmetric = TRUE,
    only.values = TRUE
  )$values
  list(
    min_eigenvalue = min(eigenvalues),
    is_psd = min(eigenvalues) >= -1e-10
  )
}

matrix_to_long <- function(value, matrix_name, replicate) {
  grid <- expand.grid(
    row = rownames(value),
    col = colnames(value),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  data.table::data.table(
    replicate = replicate,
    matrix = matrix_name,
    row = grid$row,
    col = grid$col,
    value = as.numeric(value)
  )
}

load_completed_replicate <- function(out_dir, replicate, summary_path = NULL) {
  rep_dir <- file.path(out_dir, sprintf("rep%03d", replicate))
  rds_path <- file.path(rep_dir, "background_parameters.rds")
  qc_path <- file.path(rep_dir, "qc.tsv")
  if (!file.exists(rds_path) || !file.exists(qc_path)) return(NULL)
  object <- readRDS(rds_path)
  if (is.null(object$bgparas) || is.null(object$C) || is.null(object$Omega0)) {
    return(NULL)
  }
  if (!is.null(summary_path)) {
    current_summary <- normalizePath(summary_path, mustWork = TRUE)
    current_md5 <- unname(tools::md5sum(summary_path))
    if (!identical(object$input_summary, current_summary) ||
        !identical(unname(object$input_md5), current_md5) ||
        length(current_md5) != 1L || is.na(current_md5)) return(NULL)
  }
  cross_population_intercept <- object$cross_population_intercept
  if (!(is.double(cross_population_intercept) ||
        is.integer(cross_population_intercept)) ||
      length(cross_population_intercept) != 1L ||
      !is.finite(cross_population_intercept) ||
      cross_population_intercept != 0 ||
      any(!is.finite(object$C)) ||
      any(object$C[TRAITS[1], TRAITS[2:3]] != 0) ||
      any(object$C[TRAITS[2:3], TRAITS[1]] != 0)) return(NULL)
  qc <- data.table::fread(qc_path)
  if (nrow(qc) != 1L || !is.finite(qc$m_background[[1L]])) return(NULL)
  list(
    pairs = data.table::as.data.table(object$bgparas),
    matrices = data.table::rbindlist(list(
      matrix_to_long(as.matrix(object$C), "C", replicate),
      matrix_to_long(as.matrix(object$Omega0), "Omega0", replicate)
    )),
    qc = qc,
    m_background = as.integer(qc$m_background[[1L]]),
    master_seed = object$master_seed,
    replicate_seed = object$replicate_seed
  )
}

write_matrix <- function(value, path) {
  data.table::fwrite(
    data.table::as.data.table(value, keep.rownames = "trait"),
    path,
    sep = "\t"
  )
}

run_replicate <- function(path, replicate, out_dir, scenario) {
  message(sprintf("[%03d] %s", replicate, basename(path)))
  original <- read_sumstats(path)
  m_background <- sum(original$HAS_BACKGROUND == 1L)
  dat <- remove_ambiguous(original)
  if (nrow(dat) < 1000L) stop("Too few SNPs after allele filtering")

  temporary <- file.path(out_dir, sprintf(".tmp_rep%03d", replicate))
  unlink(temporary, recursive = TRUE, force = TRUE)
  dir.create(temporary, recursive = TRUE)
  on.exit(unlink(temporary, recursive = TRUE, force = TRUE), add = TRUE)

  z_values <- write_xmr_inputs(dat, temporary)
  replicate_seed <- derive_replicate_seed(scenario$master_seed, replicate)
  set.seed(replicate_seed)
  bgparas <- estimate_bgparas_crosspop_zero(dat, z_values)
  bgparas <- standardize_bgparas(bgparas, replicate)
  bgparas <- add_stage_counts(bgparas, z_values)

  qc_input <- as.data.frame(
    bgparas[, setdiff(names(bgparas), c("replicate", "n_stage1", "n_stage2")),
              with = FALSE]
  )
  names(qc_input)[match(
    c("omega11", "omega12", "omega22"), names(qc_input)
  )] <- c("Omega11", "Omega12", "Omega22")
  matrices <- XMR::qc_bgparas(
    bgparas = qc_input,
    exposure_code = "X",
    outcome_code = "Y",
    pop_ref = "EUR",
    pop_target = "EAS",
    truncation_corr = 0.95
  )
  intercept <- matrices$C[TRAITS, TRAITS, drop = FALSE]
  stopifnot(all(intercept[1, 2:3] == 0), all(intercept[2:3, 1] == 0))
  omega <- matrices$Omega0[TRAITS, TRAITS, drop = FALSE]
  if (any(!is.finite(intercept)) || any(!is.finite(omega))) {
    stop("Non-finite C or Omega0 in replicate ", replicate)
  }

  rep_dir <- file.path(out_dir, sprintf("rep%03d", replicate))
  dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(
    bgparas, file.path(rep_dir, "background_pairs.tsv"), sep = "\t"
  )
  write_matrix(intercept, file.path(rep_dir, "C.tsv"))
  write_matrix(omega, file.path(rep_dir, "Omega0.tsv"))
  saveRDS(
    list(
      scenario = scenario,
      master_seed = scenario$master_seed,
      replicate_seed = replicate_seed,
      input_summary = normalizePath(path, mustWork = TRUE),
      input_md5 = unname(tools::md5sum(path)),
      cross_population_intercept = 0,
      bgparas = bgparas,
      C = intercept,
      Omega0 = omega
    ),
    file.path(rep_dir, "background_parameters.rds")
  )

  c_diagnostics <- matrix_diagnostics(intercept)
  omega_diagnostics <- matrix_diagnostics(omega)
  qc <- data.table::data.table(
    replicate = replicate,
    beta = scenario$beta,
    rho = scenario$rho,
    n_aux = scenario$n_aux,
    input_file = normalizePath(path),
    n_snps_input = nrow(original),
    n_ambiguous_removed = nrow(original) - nrow(dat),
    n_snps_ldsc = nrow(dat),
    m_background = m_background,
    master_seed = scenario$master_seed,
    replicate_seed = replicate_seed,
    z2_threshold = Z2_THRESHOLD,
    C_min_eigenvalue = c_diagnostics$min_eigenvalue,
    C_is_psd = c_diagnostics$is_psd,
    Omega0_min_eigenvalue = omega_diagnostics$min_eigenvalue,
    Omega0_is_psd = omega_diagnostics$is_psd,
    Omega0_negative_diagonal = any(diag(omega) < 0)
  )
  data.table::fwrite(qc, file.path(rep_dir, "qc.tsv"), sep = "\t")

  list(
    pairs = bgparas,
    matrices = data.table::rbindlist(list(
      matrix_to_long(intercept, "C", replicate),
      matrix_to_long(omega, "Omega0", replicate)
    )),
    qc = qc,
    m_background = m_background
  )
}

derive_evaluation_estimates <- function(matrices) {
  c_estimates <- data.table::copy(matrices[matrix == "C"])
  omega <- data.table::copy(matrices[matrix == "Omega0"])

  row_variances <- omega[row == col, .(
    replicate, row, variance_row = value
  )]
  col_variances <- omega[row == col, .(
    replicate, col = row, variance_col = value
  )]
  genetic_correlation <- merge(
    omega, row_variances,
    by = c("replicate", "row"), all.x = TRUE, sort = FALSE
  )
  genetic_correlation <- merge(
    genetic_correlation, col_variances,
    by = c("replicate", "col"), all.x = TRUE, sort = FALSE
  )
  genetic_correlation[, value := ifelse(
    variance_row > 0 & variance_col > 0,
    value / sqrt(variance_row * variance_col),
    NA_real_
  )]
  genetic_correlation[, matrix := "GeneticCorrelation"]

  keep <- c("replicate", "matrix", "row", "col", "value")
  data.table::rbindlist(list(
    c_estimates[, ..keep],
    omega[, ..keep],
    genetic_correlation[, ..keep]
  ))
}

make_truth_table <- function(
    m_background, rho, h2_eur_x, h2_eas_x, h2_eas_y, replicate) {
  r_omega_true <- make_r_omega_true(rho)
  result <- data.table::rbindlist(list(
    matrix_to_long(C_TRUE, "C", replicate),
    matrix_to_long(
      make_omega_true(
        m_background, rho, h2_eur_x, h2_eas_x, h2_eas_y
      ),
      "Omega0",
      replicate
    ),
    matrix_to_long(r_omega_true, "GeneticCorrelation", replicate)
  ))
  result[, `:=`(
    truth = value,
    value = NULL,
    m_background = m_background
  )]
  result
}

upper_triangle <- function(dat) {
  dat[match(row, TRAITS) <= match(col, TRAITS)]
}

evaluate_parameters <- function(estimates, truth) {
  observed <- upper_triangle(data.table::copy(estimates))
  expected <- upper_triangle(data.table::copy(truth))
  observed <- observed[matrix != "GeneticCorrelation" | row != col]
  expected <- expected[matrix != "GeneticCorrelation" | row != col]
  merged <- merge(
    observed,
    expected[, .(replicate, matrix, row, col, truth, m_background)],
    by = c("replicate", "matrix", "row", "col"),
    all.x = TRUE
  )
  merged[, element := paste(row, col, sep = "/")]

  accuracy <- merged[, .(
    truth = mean(truth),
    truth_min = min(truth),
    truth_max = max(truth),
    n_valid = sum(is.finite(value)),
    mean_estimate = mean(value, na.rm = TRUE),
    sd_estimate = stats::sd(value, na.rm = TRUE),
    median_estimate = stats::median(value, na.rm = TRUE),
    q025 = as.numeric(stats::quantile(value, 0.025, na.rm = TRUE)),
    q975 = as.numeric(stats::quantile(value, 0.975, na.rm = TRUE)),
    bias = mean(value - truth, na.rm = TRUE),
    relative_bias = if (all(abs(truth) > sqrt(.Machine$double.eps))) {
      mean((value - truth) / truth, na.rm = TRUE)
    } else {
      NA_real_
    },
    rmse = sqrt(mean((value - truth)^2, na.rm = TRUE))
  ), by = .(matrix, row, col, element)]

  merged[, weight := ifelse(row == col, 1, 2)]
  errors <- merged[, .(
    element_rmse = sqrt(mean((value - truth)^2, na.rm = TRUE)),
    frobenius_error = sqrt(sum(weight * (value - truth)^2, na.rm = TRUE)),
    relative_frobenius_error = sqrt(
      sum(weight * (value - truth)^2, na.rm = TRUE) /
        sum(weight * truth^2, na.rm = TRUE)
    )
  ), by = .(replicate, matrix)]
  list(estimates = merged, accuracy = accuracy, errors = errors)
}

draw_boxplots <- function(evaluation, out_dir) {
  matrix_order <- c("C", "Omega0", "GeneticCorrelation")
  plot_data <- data.table::copy(evaluation$estimates)
  plot_data[, matrix := factor(matrix, levels = matrix_order)]
  plot_data[, ratio := ifelse(
    abs(truth) > sqrt(.Machine$double.eps), value / truth, NA_real_
  )]
  truth_points <- plot_data[, .(truth = mean(truth)), by = .(matrix, element)]

  estimate_plot <- ggplot2::ggplot(
    plot_data, ggplot2::aes(x = element, y = value)
  ) +
    ggplot2::geom_boxplot(
      width = 0.62, outlier.shape = NA, fill = "#8EC5D6", na.rm = TRUE
    ) +
    ggplot2::geom_jitter(
      width = 0.12, alpha = 0.35, size = 1, na.rm = TRUE
    ) +
    ggplot2::geom_point(
      data = truth_points,
      ggplot2::aes(x = element, y = truth),
      inherit.aes = FALSE,
      color = "#C7352C", shape = 18, size = 3
    ) +
    ggplot2::facet_wrap(~matrix, scales = "free_y", ncol = 1L) +
    ggplot2::labs(
      title = "XMR background-parameter estimates",
      subtitle = "Red diamonds: mean truth; accuracy and ratios use replicate-specific truth",
      x = NULL, y = "Estimate"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
      panel.grid.minor = ggplot2::element_blank()
    )

  ratio_plot <- ggplot2::ggplot(
    plot_data, ggplot2::aes(x = element, y = ratio)
  ) +
    ggplot2::geom_hline(
      yintercept = 1, color = "#C7352C", linetype = "dashed"
    ) +
    ggplot2::geom_boxplot(
      width = 0.62, outlier.shape = NA, fill = "#B8C5A5", na.rm = TRUE
    ) +
    ggplot2::geom_jitter(
      width = 0.12, alpha = 0.35, size = 1, na.rm = TRUE
    ) +
    ggplot2::facet_wrap(~matrix, scales = "free_y", ncol = 1L) +
    ggplot2::labs(
      title = "Estimates relative to simulation truth",
      x = NULL, y = "Estimate / truth"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
      panel.grid.minor = ggplot2::element_blank()
    )

  ggplot2::ggsave(
    file.path(out_dir, "background_parameter_boxplots.png"),
    estimate_plot, width = 10.5, height = 12, dpi = 300
  )
  ggplot2::ggsave(
    file.path(out_dir, "background_parameter_boxplots.pdf"),
    estimate_plot, width = 10.5, height = 12
  )
  ggplot2::ggsave(
    file.path(out_dir, "background_parameter_ratio_boxplots.png"),
    ratio_plot, width = 10.5, height = 12, dpi = 300
  )
  ggplot2::ggsave(
    file.path(out_dir, "background_parameter_ratio_boxplots.pdf"),
    ratio_plot, width = 10.5, height = 12
  )
}

run_scenario <- function(
    root, beta, rho, n_aux, n_target, n_anchors,
    h2_eur_y, h2_eur_x, h2_eas_x, h2_eas_y,
    master_seed, expected_replicates, overwrite) {
  paths <- scenario_paths(
    root, beta, rho, n_aux, n_target, n_anchors,
    h2_eur_y, h2_eur_x, h2_eas_x, h2_eas_y
  )
  scenario <- list(
    setting_id = paths$setting_id,
    beta = beta,
    rho = rho,
    n_aux = n_aux,
    n_target = n_target,
    n_anchors = n_anchors,
    h2_eur_y = h2_eur_y,
    h2_eur_x = h2_eur_x,
    h2_eas_x = h2_eas_x,
    h2_eas_y = h2_eas_y,
    master_seed = master_seed
  )
  message(
    sprintf(
      "\nScenario: beta=%g, rho=%g, auxiliary N=%d, target N=%d",
      beta, rho, n_aux, n_target
    )
  )
  files <- Sys.glob(paths$input_glob)
  if (length(files) != expected_replicates) {
    stop(
      "Expected ", expected_replicates, " Step 5 files, found ",
      length(files), ": ", paths$input_glob
    )
  }
  replicates <- vapply(files, parse_replicate, integer(1))
  order_index <- order(replicates)
  files <- files[order_index]
  replicates <- replicates[order_index]
  if (anyDuplicated(replicates)) stop("Duplicated replicate numbers")

  dir.create(paths$out_dir, recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(
    data.table::data.table(
      beta = beta,
      rho = rho,
      n_aux = n_aux,
      n_target = n_target,
      n_anchors = n_anchors,
      setting_id = paths$setting_id,
      input_glob = paths$input_glob,
      output_directory = paths$out_dir,
      expected_replicates = expected_replicates,
      z2_threshold = Z2_THRESHOLD,
      h2_eur_y = h2_eur_y,
      h2_eur_x = h2_eur_x,
      h2_eas_x = h2_eas_x,
      h2_eas_y = h2_eas_y,
      master_seed = master_seed
    ),
    file.path(paths$out_dir, "scenario.tsv"),
    sep = "\t"
  )
  results <- vector("list", length(files))
  for (i in seq_along(files)) {
    completed <- if (overwrite) NULL else {
      load_completed_replicate(paths$out_dir, replicates[[i]], files[[i]])
    }
    if (!is.null(completed) &&
        (!"input_file" %in% names(completed$qc) ||
         !identical(completed$qc$input_file[[1L]], normalizePath(files[[i]])))) {
      stop("Saved background parameters belong to different summary data. ",
           "Use a new --root for the new architecture, or rerun with --overwrite.")
    }
    if (!is.null(completed) &&
        (is.null(completed$master_seed) ||
         !identical(as.integer(completed$master_seed), as.integer(master_seed)))) {
      stop("Saved background parameters use a different master seed. ",
           "Use a new --root or rerun with --overwrite.")
    }
    if (is.null(completed)) {
      results[[i]] <- run_replicate(
        files[[i]], replicates[[i]], paths$out_dir, scenario
      )
    } else {
      message(sprintf("[%03d] using existing background checkpoint", replicates[[i]]))
      results[[i]] <- completed
    }
  }

  m_background_values <- vapply(
    results, `[[`, integer(1), "m_background"
  )
  if (any(m_background_values <= 0L)) stop("M_background must be positive in every replicate")
  all_pairs <- data.table::rbindlist(lapply(results, `[[`, "pairs"))
  all_matrices <- data.table::rbindlist(lapply(results, `[[`, "matrices"))
  all_qc <- data.table::rbindlist(lapply(results, `[[`, "qc"))

  data.table::fwrite(
    all_pairs,
    file.path(paths$out_dir, "all_replicates_bgparas_pairs.tsv.gz"),
    sep = "\t"
  )
  data.table::fwrite(
    all_matrices,
    file.path(paths$out_dir, "all_replicates_matrices_long.tsv.gz"),
    sep = "\t"
  )
  data.table::fwrite(
    all_qc,
    file.path(paths$out_dir, "all_replicates_qc.tsv"),
    sep = "\t"
  )

  estimates <- derive_evaluation_estimates(all_matrices)
  truth <- data.table::rbindlist(lapply(seq_along(replicates), function(i) {
    make_truth_table(
      m_background_values[[i]], rho, h2_eur_x, h2_eas_x, h2_eas_y, replicates[[i]]
    )
  }))
  evaluation <- evaluate_parameters(estimates, truth)
  data.table::fwrite(
    estimates,
    file.path(paths$out_dir, "evaluation_estimates_long.tsv.gz"),
    sep = "\t"
  )
  data.table::fwrite(
    truth, file.path(paths$out_dir, "truth_matrices.tsv"), sep = "\t"
  )
  data.table::fwrite(
    evaluation$accuracy,
    file.path(paths$out_dir, "parameter_accuracy.tsv"),
    sep = "\t"
  )
  data.table::fwrite(
    evaluation$errors,
    file.path(paths$out_dir, "replicate_matrix_errors.tsv"),
    sep = "\t"
  )
  draw_boxplots(evaluation, paths$out_dir)

  cat("\nM_background range = ", min(m_background_values), " to ",
      max(m_background_values), "; each replicate uses its own Omega0 truth.\n", sep = "")
  cat("Results: ", paths$out_dir, "\n", sep = "")
  print(evaluation$accuracy)

  data.table::data.table(
    beta = beta,
    rho = rho,
    n_aux = n_aux,
    n_target = n_target,
    h2_eur_y = h2_eur_y,
    h2_eur_x = h2_eur_x,
    h2_eas_x = h2_eas_x,
    h2_eas_y = h2_eas_y,
    setting_id = paths$setting_id,
    m_background_mean = mean(m_background_values),
    m_background_min = min(m_background_values),
    m_background_max = max(m_background_values),
    output_directory = paths$out_dir
  )
}

main <- function() {
  settings <- parse_cli_args()
  if (is.null(settings)) return(invisible(NULL))
  require_simulation_packages()

  scenario_index <- vector(
    "list",
    length(settings$betas) * length(settings$rhos) *
      length(settings$n_aux) * length(settings$n_target)
  )
  position <- 1L
  for (beta in settings$betas) {
    for (rho in settings$rhos) {
      for (n_aux in settings$n_aux) {
        for (n_target in settings$n_target) {
          scenario_index[[position]] <- run_scenario(
            settings$root, beta, rho, n_aux, n_target,
            settings$n_anchors,
            settings$h2_eur_y, settings$h2_eur_x,
            settings$h2_eas_x, settings$h2_eas_y,
            settings$master_seed,
            settings$expected_replicates, settings$overwrite
          )
          position <- position + 1L
        }
      }
    }
  }
  index <- data.table::rbindlist(scenario_index)
  index_path <- file.path(
    settings$root, "xmr_background", "scenario_index.tsv"
  )
  dir.create(dirname(index_path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(index, index_path, sep = "\t")
  cat("\nCompleted ", nrow(index), " scenario(s).\n", sep = "")
  cat("Scenario index: ", index_path, "\n", sep = "")
}

if (sys.nframe() == 0L) main()
