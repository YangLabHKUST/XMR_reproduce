options(stringsAsFactors = FALSE, warn = 1)

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x)) y else x
}

get_arg <- function(args, key, default = NULL) {
  hit <- grep(paste0("^", key, "="), args)
  if (length(hit)) return(sub(paste0("^", key, "="), "", args[hit[1L]]))
  hit <- which(args == key)
  if (length(hit) && hit[1L] < length(args)) return(args[hit[1L] + 1L])
  default
}

has_flag <- function(args, key) {
  any(args == key)
}

derive_simulation_seed <- function(
    master_seed, component_id, replicate = 0L, analysis_index = 0L) {
  values <- c(master_seed, component_id, replicate, analysis_index)
  if (any(!is.finite(values)) || any(values < 0) || any(values != floor(values))) {
    stop("Seed inputs must be non-negative integers", call. = FALSE)
  }
  modulus <- 2147483646
  mixed <- (as.double(master_seed) +
    1000003 * as.double(component_id) +
    9176 * as.double(replicate) +
    6113 * as.double(analysis_index)) %% modulus
  as.integer(mixed + 1)
}

env_or <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

threshold_label <- function(value) {
  sub("e([+-])0+", "e\\1", formatC(value, format = "e", digits = 0L))
}

DEFAULT_WORK_ROOT <- paste0(
  "simulations/data/",
  "reduced_confounding"
)
DEFAULT_DATA_ROOT <- "simulations/data"

parse_thresholds <- function(text) {
  values <- suppressWarnings(as.numeric(strsplit(text, ",", fixed = TRUE)[[1L]]))
  if (!length(values) || any(!is.finite(values)) || any(values <= 0 | values >= 1)) {
    stop("--thresholds must be comma-separated values in (0, 1)", call. = FALSE)
  }
  unique(values)
}

parse_method_config <- function(args, method, script_dir) {
  root <- get_arg(
    args,
    "--root",
    env_or("SIM_ROOT", DEFAULT_WORK_ROOT)
  )
  data_root <- get_arg(
    args, "--data-root", env_or("SIM_DATA_ROOT", DEFAULT_DATA_ROOT)
  )
  setting_id <- get_arg(args, "--setting-id", env_or("SIM_SETTING_ID", ""))
  if (!nzchar(setting_id)) {
    stop("--setting-id is required", call. = FALSE)
  }
  output_root <- get_arg(
    args,
    "--output-root",
    env_or("SIM_METHOD_OUTPUT_ROOT", file.path(root, "method_results"))
  )
  replicate_text <- get_arg(
    args,
    "--replicate",
    env_or("SIM_REPLICATE", "")
  )
  replicate <- if (nzchar(replicate_text)) {
    suppressWarnings(as.integer(replicate_text))
  } else {
    NA_integer_
  }
  if (nzchar(replicate_text) && (is.na(replicate) || replicate < 1L)) {
    stop("--replicate must be a positive integer", call. = FALSE)
  }
  default_output_dir <- file.path(output_root, setting_id, method)
  if (!is.na(replicate)) {
    default_output_dir <- file.path(
      default_output_dir, sprintf("rep%03d", replicate)
    )
  }

  input_glob <- get_arg(args, "--input-glob", env_or("SIM_INPUT_GLOB", ""))
  if (!nzchar(input_glob)) {
    stop("--input-glob is required", call. = FALSE)
  }
  background_dir <- get_arg(
    args,
    "--background-dir",
    env_or("SIM_BACKGROUND_DIR", file.path(root, "xmr_background", setting_id))
  )
  common_snps <- get_arg(
    args, "--common-snps", env_or("SIM_COMMON_SNPS", "")
  )
  if (!nzchar(common_snps)) {
    matches <- sort(Sys.glob(input_glob))
    if (length(matches)) common_snps <- matches[[1L]]
  }

  cfg <- list(
    method = method,
    setting_id = setting_id,
    root = normalizePath(root, mustWork = FALSE),
    data_root = normalizePath(data_root, mustWork = FALSE),
    input_glob = input_glob,
    background_dir = background_dir,
    aux_background_dir = get_arg(
      args,
      "--aux-background-dir",
      file.path(root, "mrbee_aux_background", setting_id)
    ),
    output_dir = get_arg(args, "--output-dir", default_output_dir),
    common_snps = common_snps,
    eur_schema = get_arg(
      args,
      "--eur-schema",
      file.path(data_root, "vilma_ld/wb_hapmap/ld_manifest.txt")
    ),
    eas_schema = get_arg(
      args,
      "--eas-schema",
      file.path(data_root, "vilma_ld/eas_hapmap/ld_manifest.txt")
    ),
    clumper = get_arg(
      args,
      "--clumper",
      file.path(script_dir, "vilma_ld_clump.py")
    ),
    python = get_arg(args, "--python", env_or("SIM_PYTHON", "python3")),
    expected_replicates = suppressWarnings(as.integer(get_arg(
      args, "--expected-replicates", "50"
    ))),
    expected_snps = suppressWarnings(as.integer(get_arg(
      args, "--expected-snps", "0"
    ))),
    replicate = replicate,
    thresholds = parse_thresholds(get_arg(args, "--thresholds", "5e-5,5e-8")),
    clump_kb = suppressWarnings(as.integer(get_arg(args, "--clump-kb", "10000"))),
    clump_r2 = suppressWarnings(as.numeric(get_arg(args, "--clump-r2", "0.001"))),
    beta_true = suppressWarnings(as.numeric(get_arg(args, "--beta-true", "0"))),
    rho = suppressWarnings(as.numeric(get_arg(args, "--rho", "NA"))),
    n_aux = suppressWarnings(as.integer(get_arg(args, "--n-aux", "NA"))),
    n_target = suppressWarnings(as.integer(get_arg(args, "--n-target", "NA"))),
    aux_z2_threshold = suppressWarnings(as.numeric(get_arg(
      args, "--aux-z2-threshold", "30"
    ))),
    seed = suppressWarnings(as.integer(get_arg(args, "--seed", "414"))),
    overwrite = has_flag(args, "--overwrite") ||
      tolower(env_or("SIM_OVERWRITE", "false")) %in% c("1", "true", "yes")
  )
  validate_method_config(cfg)
  cfg
}

validate_method_config <- function(cfg) {
  if (is.na(cfg$seed) || cfg$seed < 0L || cfg$seed > 2147483646L) {
    stop("--seed must be an integer between 0 and 2147483646", call. = FALSE)
  }
  if (is.na(cfg$expected_replicates) || cfg$expected_replicates < 0L) {
    stop("--expected-replicates must be non-negative", call. = FALSE)
  }
  if (is.na(cfg$expected_snps) || cfg$expected_snps < 0L) {
    stop("--expected-snps must be non-negative", call. = FALSE)
  }
  if (is.na(cfg$clump_kb) || cfg$clump_kb <= 0L) {
    stop("--clump-kb must be positive", call. = FALSE)
  }
  if (!is.finite(cfg$clump_r2) || cfg$clump_r2 <= 0 || cfg$clump_r2 > 1) {
    stop("--clump-r2 must be in (0, 1]", call. = FALSE)
  }
  if (!is.finite(cfg$beta_true)) stop("--beta-true must be finite", call. = FALSE)
  if (!is.finite(cfg$rho) || abs(cfg$rho) >= 1) {
    stop("--rho must be finite and in (-1, 1)", call. = FALSE)
  }
  if (is.na(cfg$n_aux) || cfg$n_aux <= 0L) {
    stop("--n-aux must be positive", call. = FALSE)
  }
  if (is.na(cfg$n_target) || cfg$n_target <= 0L) {
    stop("--n-target must be positive", call. = FALSE)
  }
  if (!is.finite(cfg$aux_z2_threshold) || cfg$aux_z2_threshold <= 0) {
    stop("--aux-z2-threshold must be positive", call. = FALSE)
  }
  if (!nzchar(cfg$common_snps)) {
    stop("No summary file is available for --common-snps", call. = FALSE)
  }
  invisible(cfg)
}

require_method_packages <- function(packages) {
  missing <- packages[!vapply(
    packages, requireNamespace, logical(1), quietly = TRUE
  )]
  if (length(missing)) {
    stop("Missing R package(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

parse_replicate <- function(path) {
  name <- basename(path)
  match <- regexec("_rep([0-9]+)\\.tsv\\.gz$", name)
  fields <- regmatches(name, match)[[1L]]
  if (length(fields) != 2L) {
    stop("Cannot parse replicate from summary filename: ", path, call. = FALSE)
  }
  as.integer(fields[2L])
}

discover_summary_files <- function(cfg) {
  files <- sort(Sys.glob(cfg$input_glob))
  if (!length(files)) stop("No summary files matched: ", cfg$input_glob, call. = FALSE)
  replicates <- vapply(files, parse_replicate, integer(1))
  if (anyDuplicated(replicates)) stop("Duplicated replicate numbers", call. = FALSE)
  ordering <- order(replicates)
  files <- files[ordering]
  replicates <- replicates[ordering]

  if (!is.na(cfg$replicate)) {
    keep <- replicates == cfg$replicate
    if (!any(keep)) stop("Requested replicate was not found: ", cfg$replicate)
    files <- files[keep]
    replicates <- replicates[keep]
  } else if (cfg$expected_replicates > 0L &&
             length(files) != cfg$expected_replicates) {
    stop(
      "Expected ", cfg$expected_replicates, " summary files but found ",
      length(files), call. = FALSE
    )
  }
  data.table::data.table(replicate = replicates, path = files)
}

read_tsv <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE) &&
      !requireNamespace("R.utils", quietly = TRUE)) {
    connection <- gzfile(path, open = "rt")
    on.exit(close(connection), add = TRUE)
    return(data.table::as.data.table(utils::read.delim(
      connection,
      header = TRUE,
      sep = "\t",
      check.names = FALSE,
      stringsAsFactors = FALSE
    )))
  }
  data.table::fread(path, showProgress = FALSE)
}

read_simulation_summary <- function(path, expected_snps = 0L) {
  dat <- read_tsv(path)
  required <- c(
    "SNP", "CHR", "BP", "A1", "A2",
    "b.out.pop1", "se.out.pop1", "pval.out.pop1", "N.out.pop1",
    "b.exp.pop1", "se.exp.pop1", "pval.exp.pop1", "N.exp.pop1",
    "b.exp.pop2", "se.exp.pop2", "pval.exp.pop2", "N.exp.pop2",
    "b.out.pop2", "se.out.pop2", "pval.out.pop2", "N.out.pop2",
    "L2.pop1", "L2.pop2", "L12"
  )
  missing <- setdiff(required, names(dat))
  if (length(missing)) {
    stop(path, " is missing columns: ", paste(missing, collapse = ", "))
  }
  if (expected_snps > 0L && nrow(dat) != expected_snps) {
    stop("Expected ", expected_snps, " SNPs but found ", nrow(dat), " in ", path)
  }
  dat[, `:=`(
    SNP = as.character(SNP),
    A1 = toupper(as.character(A1)),
    A2 = toupper(as.character(A2)),
    CHR = as.character(CHR),
    BP = as.integer(BP)
  )]
  if (anyNA(dat$SNP) || any(!nzchar(dat$SNP)) || anyDuplicated(dat$SNP)) {
    stop("Invalid or duplicated SNP IDs in ", path)
  }
  numeric_columns <- setdiff(required, c("SNP", "CHR", "A1", "A2"))
  for (column in numeric_columns) {
    data.table::set(dat, j = column, value = suppressWarnings(as.numeric(dat[[column]])))
  }
  if (any(!is.finite(as.matrix(dat[, ..numeric_columns])))) {
    stop("Non-finite required numeric values in ", path)
  }
  dat
}

background_object_md5 <- function(object) {
  path <- tempfile("background_signature_", fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(object, path, version = 2)
  unname(as.character(tools::md5sum(path)))
}

background_file_identity <- function(path, required = FALSE) {
  path <- path.expand(as.character(path))
  if (length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("Background provenance requires one file path")
  }
  if (!grepl("^/", path)) path <- file.path(getwd(), path)
  normalized <- normalizePath(path, mustWork = FALSE)
  exists <- file.exists(path) && !dir.exists(path)
  if (required && !exists) stop("Missing background input: ", path)
  md5 <- if (exists) unname(as.character(tools::md5sum(path))) else NA_character_
  if (exists && (length(md5) != 1L || is.na(md5))) stop("Cannot hash background input: ", path)
  list(path = normalized, exists = exists, md5 = md5)
}

background_package_versions <- function(packages) {
  stats::setNames(lapply(packages, function(package) {
    tryCatch(as.character(utils::packageVersion(package)), error = function(e) NA_character_)
  }), packages)
}

analysis_fingerprint <- function(
    summary_path, background_paths, code_paths = character(),
    parameters = list(), packages = character()) {
  identity <- list(
    schema = "method-background-provenance-v1",
    summary = background_file_identity(summary_path),
    backgrounds = lapply(background_paths, background_file_identity),
    code = lapply(code_paths, background_file_identity),
    parameters = parameters,
    packages = background_package_versions(packages)
  )
  paste0("analysis-v1:", background_object_md5(identity))
}

normalize_background_matrix <- function(value, traits, label, path) {
  value <- as.matrix(value)
  storage.mode(value) <- "double"
  if (nrow(value) != ncol(value) || !nrow(value) %in% c(2L, 3L)) {
    stop(label, " must be a 2 x 2 or 3 x 3 matrix in ", path)
  }
  if (is.null(rownames(value)) && is.null(colnames(value))) {
    order <- if (nrow(value) == 3L) c("EUR_X", "EAS_X", "EAS_Y") else traits
    if (length(order) != nrow(value)) stop(label, " has ambiguous trait order in ", path)
    dimnames(value) <- list(order, order)
  } else if (is.null(rownames(value)) || is.null(colnames(value))) {
    stop(label, " must provide both row and column trait names in ", path)
  }
  if (anyDuplicated(rownames(value)) || anyDuplicated(colnames(value))) {
    stop(label, " contains duplicated trait names in ", path)
  }
  missing <- setdiff(traits, intersect(rownames(value), colnames(value)))
  if (length(missing)) stop(label, " lacks traits: ", paste(missing, collapse = ", "))
  value <- value[traits, traits, drop = FALSE]
  if (any(!is.finite(value))) stop(label, " contains non-finite values in ", path)
  0.5 * (value + t(value))
}

load_background_parameters <- function(
    cfg, replicate, require_crosspop_zero = TRUE, summary_path = NULL) {
  path <- file.path(cfg$background_dir, sprintf("rep%03d", replicate), "background_parameters.rds")
  if (!file.exists(path)) stop("Missing background parameter file: ", path)
  object <- readRDS(path)
  if (is.null(object$C) || is.null(object$Omega0)) stop(path, " must contain C and Omega0")
  if (!is.null(object$master_seed) && !is.null(cfg$seed) &&
      !identical(as.integer(object$master_seed), as.integer(cfg$seed))) {
    stop("Background checkpoint uses a different master seed: ", path)
  }
  for (field in c("beta", "rho", "n_aux", "n_target")) {
    saved <- object$scenario[[field]]
    requested <- if (field == "beta") cfg$beta_true else cfg[[field]]
    if (!is.null(saved) && !is.null(requested) &&
        !isTRUE(all.equal(as.numeric(saved), as.numeric(requested), check.attributes = FALSE))) {
      stop("Background checkpoint uses a different ", field, ": ", path)
    }
  }
  if (!is.null(summary_path) && !is.null(object$input_md5)) {
    input <- background_file_identity(summary_path, required = TRUE)
    if (!identical(object$input_summary, input$path) ||
        !identical(as.character(object$input_md5), input$md5)) {
      stop("Background checkpoint belongs to different summary data. Rerun step6: ", path)
    }
  }
  # MRBEE/MR-APSS consume only this within-EAS block. Cross-population
  # intercept requirements belong to callers that actually use the full matrix.
  traits <- if (require_crosspop_zero) c("EUR_X", "EAS_X", "EAS_Y") else c("EAS_X", "EAS_Y")
  C <- normalize_background_matrix(object$C, traits, "C", path)
  Omega0 <- normalize_background_matrix(object$Omega0, traits, "Omega0", path)
  if (any(diag(C) <= 0) || min(eigen(C, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
    stop("Background C is not positive definite: ", path)
  }
  if (require_crosspop_zero) {
    marker <- object$cross_population_intercept
    valid_marker <- is.numeric(marker) && length(marker) == 1L && is.finite(marker) && marker == 0
    if (!valid_marker || any(C[1, 2:3] != 0) || any(C[2:3, 1] != 0)) {
      stop("Background checkpoint predates the zero cross-population intercept fit. Rerun step6: ", path)
    }
  }
  list(C = C, Omega0 = Omega0, path = path)
}

target_background <- function(background) {
  traits <- c("EAS_X", "EAS_Y")
  list(
    C = background$C[traits, traits, drop = FALSE],
    Omega0 = background$Omega0[traits, traits, drop = FALSE]
  )
}

estimate_gc_pair <- function(z1, n1, z2, n2, l1, l2, l12, z2_threshold = 30) {
  keep <- is.finite(z1) & is.finite(z2) &
    is.finite(n1) & n1 > 0 & is.finite(n2) & n2 > 0 &
    is.finite(l1) & is.finite(l2) & is.finite(l12) &
    z1^2 < z2_threshold & z2^2 < z2_threshold
  if (sum(keep) < 1000L) stop("Too few SNPs remain for source-population LDSC")
  gwas1 <- data.frame(Z = z1, N = n1)
  gwas2 <- data.frame(Z = z2, N = n2)
  w1 <- 1 / pmax(l1, 1)
  w2 <- 1 / pmax(l2, 1)
  wx <- sqrt(w1 * w2)
  fit1 <- XMAP::estimate_gc(
    gwas1[keep, , drop = FALSE],
    gwas2[keep, , drop = FALSE],
    l1[keep], l2[keep], l12[keep],
    reg_w1 = w1[keep], reg_w2 = w2[keep], reg_wx = wx[keep],
    constrain_intercept = FALSE
  )
  fit2 <- XMAP::estimate_gc(
    gwas1, gwas2, l1, l2, l12,
    reg_w1 = w1, reg_w2 = w2, reg_wx = wx,
    constrain_intercept = TRUE,
    fit1$tau1$coefs[1L], fit1$tau2$coefs[1L], fit1$theta$coefs[1L]
  )
  c(
    omega1 = fit2$tau1$coefs[2L],
    omega12 = fit2$theta$coefs[2L],
    omega2 = fit2$tau2$coefs[2L],
    c1 = fit2$tau1$coefs[1L],
    c12 = fit2$theta$coefs[1L],
    c2 = fit2$tau2$coefs[1L]
  )
}

cap_covariance_correlation <- function(matrix, cap = 0.95, label = "matrix") {
  matrix <- 0.5 * (matrix + t(matrix))
  if (any(!is.finite(matrix)) || any(diag(matrix) <= 0)) {
    stop(label, " has non-positive or non-finite diagonal entries")
  }
  for (pair in combn(seq_len(nrow(matrix)), 2L, simplify = FALSE)) {
    i <- pair[[1L]]
    j <- pair[[2L]]
    bound <- sqrt(matrix[i, i] * matrix[j, j])
    if (abs(matrix[i, j]) > cap * bound) {
      matrix[i, j] <- matrix[j, i] <- cap * bound * sign(matrix[i, j])
    }
  }
  if (min(eigen(matrix, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
    stop(label, " is not positive definite after correlation capping")
  }
  matrix
}

estimate_source_background <- function(
    dat, z2_threshold = 30, require_positive_omega = TRUE) {
  zx <- dat$b.exp.pop1 / dat$se.exp.pop1
  zy <- dat$b.out.pop1 / dat$se.out.pop1
  nx <- dat$N.exp.pop1
  ny <- dat$N.out.pop1
  ld <- dat$L2.pop1
  fit_x <- estimate_gc_pair(zx, nx, zx, nx, ld, ld, ld, z2_threshold)
  fit_y <- estimate_gc_pair(zy, ny, zy, ny, ld, ld, ld, z2_threshold)
  fit_xy <- estimate_gc_pair(zx, nx, zy, ny, ld, ld, ld, z2_threshold)
  traits <- c("EUR_X", "EUR_Y")
  C <- matrix(
    c(fit_x[["c1"]], fit_xy[["c12"]], fit_xy[["c12"]], fit_y[["c1"]]),
    2L,
    byrow = TRUE,
    dimnames = list(traits, traits)
  )
  Omega0 <- matrix(
    c(
      fit_x[["omega1"]], fit_xy[["omega12"]],
      fit_xy[["omega12"]], fit_y[["omega1"]]
    ),
    2L,
    byrow = TRUE,
    dimnames = list(traits, traits)
  )
  C <- cap_covariance_correlation(C, label = "source C")
  if (require_positive_omega) {
    Omega0 <- cap_covariance_correlation(Omega0, label = "source Omega0")
  }
  list(C = C, Omega0 = Omega0)
}

load_or_estimate_aux_background <- function(cfg, replicate, summary_path) {
  rep_dir <- file.path(cfg$aux_background_dir, sprintf("rep%03d", as.integer(replicate)))
  path <- file.path(rep_dir, "background_parameters.rds")
  input <- background_file_identity(summary_path, required = TRUE)
  signature <- background_object_md5(list(
    schema = "aux-background-provenance-v1", setting_id = cfg$setting_id,
    replicate = as.integer(replicate), input = input,
    z2_threshold = as.numeric(cfg$aux_z2_threshold), expected_snps = cfg$expected_snps,
    estimator = lapply(list(read_tsv, read_simulation_summary, estimate_gc_pair,
                           estimate_source_background, cap_covariance_correlation),
                      function(f) list(formals = formals(f), body = body(f))),
    packages = background_package_versions("XMAP")
  ))
  if (file.exists(path) && !cfg$overwrite) {
    object <- tryCatch(readRDS(path), error = function(e) NULL)
    if (!is.null(object) && identical(object$cache_signature, signature)) {
      checked <- tryCatch(list(
        C = normalize_background_matrix(object$C, c("EUR_X", "EUR_Y"), "source C", path),
        Omega0 = normalize_background_matrix(object$Omega0, c("EUR_X", "EUR_Y"), "source Omega0", path)
      ), error = function(e) NULL)
      if (!is.null(checked) && all(vapply(checked, function(x) {
          all(diag(x) > 0) && min(eigen(x, symmetric = TRUE, only.values = TRUE)$values) > 0
        }, logical(1)))) {
        object$C <- checked$C; object$Omega0 <- checked$Omega0
        return(object)
      }
    }
    message("Re-estimating stale auxiliary background: ", path)
  }
  dat <- read_simulation_summary(summary_path, cfg$expected_snps)
  estimate <- estimate_source_background(dat, z2_threshold = cfg$aux_z2_threshold,
                                       require_positive_omega = TRUE)
  object <- list(
    setting_id = cfg$setting_id, replicate = as.integer(replicate),
    traits = c("EUR_X", "EUR_Y"), z2_threshold = cfg$aux_z2_threshold,
    input_summary = input$path, input_md5 = input$md5, cache_signature = signature,
    C = estimate$C, Omega0 = estimate$Omega0
  )
  dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp")
  saveRDS(object, temporary)
  if (!file.rename(temporary, path)) {
    unlink(temporary)
    stop("Could not finalize auxiliary background file: ", path)
  }
  object
}

make_job_id <- function(method, replicate, threshold, mode, population = "") {
  paste(
    method,
    sprintf("rep%03d", replicate),
    threshold_label(threshold),
    mode,
    population,
    sep = "__"
  )
}

run_vilma_clump <- function(jobs, cfg, label) {
  required <- c("job_id", "population", "SNP", "priority")
  missing <- setdiff(required, names(jobs))
  if (length(missing)) stop("Clump jobs lack: ", paste(missing, collapse = ", "))
  jobs <- data.table::as.data.table(jobs)[, ..required]
  jobs[, `:=`(
    job_id = as.character(job_id),
    population = as.character(population),
    SNP = as.character(SNP),
    priority = as.numeric(priority)
  )]
  jobs <- jobs[
    nzchar(job_id) & nzchar(population) & nzchar(SNP) & is.finite(priority)
  ]
  data.table::setorder(jobs, job_id, priority, SNP)
  jobs <- unique(jobs, by = c("job_id", "SNP"))
  if (!nrow(jobs)) {
    return(data.table::data.table(
      job_id = character(), population = character(), SNP = character(),
      priority = numeric(), selection_rank = integer()
    ))
  }

  dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
  job_path <- file.path(cfg$output_dir, paste0(label, "_clump_candidates.tsv.gz"))
  selected_path <- file.path(cfg$output_dir, paste0(label, "_clump_selected.tsv.gz"))
  summary_path <- file.path(cfg$output_dir, paste0(label, "_clump_summary.tsv"))
  log_path <- file.path(cfg$output_dir, paste0(label, "_clump.log"))

  refresh_stale <- isTRUE(cfg$refresh_stale_clump)
  provenance_path <- paste0(summary_path, ".provenance.rds")
  provenance <- if (refresh_stale) analysis_fingerprint(
    cfg$common_snps, c(cfg$eur_schema, cfg$eas_schema),
    code_paths = cfg$clumper,
    parameters = list(clump_kb = cfg$clump_kb, clump_r2 = cfg$clump_r2,
                      expected_snps = cfg$expected_snps)
  ) else NULL
  if (file.exists(selected_path) && file.exists(summary_path) && !cfg$overwrite) {
    if (!file.exists(job_path) && !refresh_stale) {
      stop("Existing clump output lacks its candidate file: ", selected_path)
    }
    old_jobs <- if (file.exists(job_path)) read_tsv(job_path)[, ..required] else jobs[0]
    data.table::setorder(old_jobs, job_id, priority, SNP)
    old_jobs <- unique(old_jobs, by = c("job_id", "SNP"))
    same_identifiers <- identical(
      old_jobs[, .(job_id, population, SNP)],
      jobs[, .(job_id, population, SNP)]
    )
    same_priorities <- isTRUE(all.equal(
      old_jobs$priority, jobs$priority, tolerance = 1e-14
    ))
    same_provenance <- !refresh_stale || (file.exists(provenance_path) &&
      identical(tryCatch(readRDS(provenance_path), error = function(e) NULL), provenance))
    if (same_identifiers && same_priorities && same_provenance) return(read_tsv(selected_path))
    if (!refresh_stale) {
      stop("Existing clump cache does not match current candidates; use --overwrite")
    }
    message("Reclumping stale candidate/LD cache: ", selected_path)
  }

  # Invalidate the old transaction before replacing its candidate input. A failed
  # reclump must never associate new candidates with old selected output.
  if (refresh_stale) unlink(c(provenance_path, selected_path, summary_path))
  data.table::fwrite(jobs, job_path, sep = "\t", quote = FALSE)
  command_args <- c(
    shQuote(cfg$clumper),
    "--jobs", shQuote(job_path),
    "--output", shQuote(selected_path),
    "--summary", shQuote(summary_path),
    "--common-snps", shQuote(cfg$common_snps),
    "--eur-schema", shQuote(cfg$eur_schema),
    "--eas-schema", shQuote(cfg$eas_schema),
    "--clump-kb", as.character(cfg$clump_kb),
    "--clump-r2", format(cfg$clump_r2, scientific = FALSE),
    "--expected-snps", as.character(cfg$expected_snps),
    "--overwrite"
  )
  status <- system2(
    cfg$python,
    args = command_args,
    stdout = log_path,
    stderr = log_path
  )
  if (!identical(status, 0L) || !file.exists(selected_path)) {
    tail_log <- if (file.exists(log_path)) {
      paste(tail(readLines(log_path, warn = FALSE), 20L), collapse = "\n")
    } else {
      "No clump log was created."
    }
    stop("ViLMA clumping failed:\n", tail_log, call. = FALSE)
  }
  selected <- read_tsv(selected_path)
  selected_required <- c("job_id", "population", "SNP", "priority", "selection_rank")
  if (!all(selected_required %in% names(selected))) stop("Invalid clump output: ", selected_path)
  if (refresh_stale) saveRDS(provenance, provenance_path)
  selected
}

selected_snps <- function(selected, job_id) {
  as.character(selected$SNP[selected$job_id == job_id])
}

scalar_or_na <- function(x) {
  if (is.null(x) || !length(x)) return(NA_real_)
  value <- suppressWarnings(as.numeric(x[[1L]]))
  if (is.finite(value)) value else NA_real_
}

p_from_beta_se <- function(beta, se) {
  if (!is.finite(beta) || !is.finite(se) || se <= 0) return(NA_real_)
  2 * stats::pnorm(-abs(beta / se))
}

report_method_result <- function(result, run_index, total_runs) {
  if (nrow(result) != 1L) {
    stop("Progress reporting requires exactly one result row")
  }
  format_value <- function(value, digits = 5L) {
    value <- suppressWarnings(as.numeric(value[[1L]]))
    if (is.finite(value)) format(value, digits = digits) else "NA"
  }
  row <- result[1L]
  text <- sprintf(
    "[%d/%d] %s rep=%03d P=%s IV=%s status=%s beta=%s SE=%s p=%s",
    as.integer(run_index), as.integer(total_runs), as.character(row$method),
    as.integer(row$replicate), threshold_label(row$threshold),
    ifelse(is.na(row$nIV), "NA", as.character(row$nIV)),
    as.character(row$status), format_value(row$beta),
    format_value(row$se), format_value(row$p)
  )
  if ("nvalid" %in% names(row)) {
    text <- paste0(text, " nvalid=", format_value(row$nvalid))
  }
  origin_columns <- c(
    "n_foreground_IV", "n_shared_IV", "n_target_only_IV", "n_aux_only_IV",
    "n_background_only_IV", "n_null_effect_IV"
  )
  if (all(origin_columns %in% names(row))) {
    text <- paste0(text, sprintf(
      " foreground=%s (shared=%s target_only=%s aux_only=%s) background_only=%s null_effect=%s",
      format_value(row$n_foreground_IV), format_value(row$n_shared_IV),
      format_value(row$n_target_only_IV), format_value(row$n_aux_only_IV),
      format_value(row$n_background_only_IV), format_value(row$n_null_effect_IV)
    ))
  }
  if (!is.na(row$error) && nzchar(as.character(row$error))) {
    text <- paste0(text, " | ", as.character(row$error))
  }
  message(text)
  flush.console()
}

simulation_result <- function(
    cfg, replicate, method, threshold, selection_mode, population,
    beta = NA_real_, se = NA_real_, p = NA_real_, n_iv = NA_integer_,
    status = "ok", error = NA_character_) {
  if (!is.finite(p)) p <- p_from_beta_se(beta, se)
  ci_lower <- if (is.finite(beta) && is.finite(se) && se > 0) {
    beta - stats::qnorm(0.975) * se
  } else {
    NA_real_
  }
  ci_upper <- if (is.finite(beta) && is.finite(se) && se > 0) {
    beta + stats::qnorm(0.975) * se
  } else {
    NA_real_
  }
  data.table::data.table(
    setting_id = cfg$setting_id,
    replicate = as.integer(replicate),
    method = method,
    population = population,
    threshold = as.numeric(threshold),
    selection_mode = selection_mode,
    beta_true = cfg$beta_true,
    rho = cfg$rho,
    n_aux = cfg$n_aux,
    n_target = cfg$n_target,
    master_seed = cfg$seed,
    beta = as.numeric(beta),
    se = as.numeric(se),
    p = as.numeric(p),
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    nIV = as.integer(n_iv),
    status = status,
    error = error
  )
}

prepare_result_file <- function(cfg, filename = "results.tsv") {
  dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(cfg$output_dir, filename)
  if (file.exists(path) && cfg$overwrite) unlink(path)
  path
}

result_key <- function(replicate, method, threshold, selection_mode, population) {
  paste(
    as.integer(replicate), method, format(as.numeric(threshold), digits = 17),
    selection_mode, population, sep = "\r"
  )
}

completed_result_keys <- function(path, master_seed, analysis_signatures = NULL) {
  if (!file.exists(path)) return(character())
  existing <- read_tsv(path)
  required <- c(
    "replicate", "method", "threshold", "selection_mode", "population"
  )
  missing <- setdiff(required, names(existing))
  if (length(missing)) {
    stop("Existing result file lacks: ", paste(missing, collapse = ", "))
  }
  if (!"master_seed" %in% names(existing)) {
    stop(
      "Existing result file predates the master-seed metadata; ",
      "use a new output directory or --overwrite.", call. = FALSE
    )
  }
  saved_seeds <- unique(existing$master_seed[!is.na(existing$master_seed)])
  if (length(saved_seeds) != 1L ||
      !identical(as.integer(saved_seeds), as.integer(master_seed))) {
    stop(
      "Existing result file uses a different master seed; ",
      "use a new output directory or --overwrite.", call. = FALSE
    )
  }
  if (!is.null(analysis_signatures) &&
      !all(c("status", "beta", "se") %in% names(existing))) return(character())
  if ("status" %in% names(existing)) {
    complete <- !is.na(existing$status) & existing$status == "ok"
    if (all(c("beta", "se") %in% names(existing))) {
      complete <- complete & is.finite(existing$beta) &
        is.finite(existing$se) & existing$se > 0
    }
    existing <- existing[complete]
  }
  if (!is.null(analysis_signatures)) {
    if (!"analysis_signature" %in% names(existing)) return(character())
    if (is.null(names(analysis_signatures)) || anyDuplicated(names(analysis_signatures))) {
      stop("Analysis signatures must be named by replicate")
    }
    expected <- unname(analysis_signatures[as.character(existing$replicate)])
    keep <- !is.na(expected) & nzchar(expected) &
      !is.na(existing$analysis_signature) & existing$analysis_signature == expected
    existing <- existing[keep]
  }
  mapply(
    result_key,
    existing$replicate,
    existing$method,
    existing$threshold,
    existing$selection_mode,
    existing$population,
    USE.NAMES = FALSE
  )
}

write_checkpoint <- function(results, path, provenance_column = NULL) {
  if (!length(results)) {
    if (file.exists(path)) return(invisible(read_tsv(path)))
    return(invisible(data.table::data.table()))
  }
  incoming <- data.table::rbindlist(results, fill = TRUE, use.names = TRUE)
  if (!is.null(provenance_column) &&
      (length(provenance_column) != 1L || !provenance_column %in% names(incoming))) {
    stop("Incoming results must contain the requested provenance column")
  }
  key <- c(
    "setting_id", "replicate", "method", "threshold",
    "selection_mode", "population"
  )
  if (!nrow(incoming)) {
    if (file.exists(path)) return(invisible(read_tsv(path)))
    return(invisible(incoming))
  }
  valid_estimate <- is.finite(incoming$beta) & is.finite(incoming$se) &
    incoming$se > 0
  invalid_success <- !is.na(incoming$status) & incoming$status == "ok" &
    !valid_estimate
  incoming[invalid_success, `:=`(
    status = "failed",
    error = "Method returned status=ok without a finite beta and positive SE"
  )]
  incoming[valid_estimate, `:=`(
    p = ifelse(
      is.finite(p), p,
      2 * stats::pnorm(-abs(beta / se))
    ),
    ci_lower = beta - stats::qnorm(0.975) * se,
    ci_upper = beta + stats::qnorm(0.975) * se
  )]
  incoming <- unique(incoming, by = key)

  if (!file.exists(path)) {
    data.table::fwrite(
      incoming, path, sep = "\t", quote = FALSE, na = "NA"
    )
    return(invisible(incoming))
  }

  existing <- read_tsv(path)
  schema_migrated <- FALSE
  if (!is.null(provenance_column)) {
    if (length(provenance_column) != 1L || !provenance_column %in% names(incoming)) {
      stop("Incoming results must contain the requested provenance column")
    }
    if (!provenance_column %in% names(existing)) {
      data.table::set(existing, j = provenance_column, value = NA_character_)
      schema_migrated <- TRUE
    }
  }
  missing_existing <- setdiff(names(incoming), names(existing))
  missing_incoming <- setdiff(names(existing), names(incoming))
  if (length(missing_existing) || length(missing_incoming)) {
    stop(
      "Result schema changed for ", path,
      "; use --overwrite. Existing-only: ",
      paste(missing_incoming, collapse = ","),
      "; incoming-only: ", paste(missing_existing, collapse = ",")
    )
  }
  existing_keys <- do.call(paste, c(existing[, ..key], sep = "\r"))
  incoming_keys <- do.call(paste, c(incoming[, ..key], sep = "\r"))

  # Retry failed rows and replace successes produced under different provenance.
  # Matching successful rows remain unchanged.
  existing_complete <- !is.na(existing$status) & existing$status == "ok"
  if (all(c("beta", "se") %in% names(existing))) {
    existing_complete <- existing_complete & is.finite(existing$beta) &
      is.finite(existing$se) & existing$se > 0
  }
  incoming_match <- match(existing_keys, incoming_keys)
  provenance_changed <- rep(FALSE, nrow(existing))
  if (!is.null(provenance_column)) {
    expected <- incoming[[provenance_column]][incoming_match]
    saved <- existing[[provenance_column]]
    provenance_changed <- is.na(saved) | is.na(expected) | saved != expected
  }
  retry_rows <- which(!is.na(incoming_match) & (!existing_complete | provenance_changed))
  if (length(retry_rows)) {
    existing <- existing[-retry_rows]
    existing_keys <- existing_keys[-retry_rows]
  }

  new_rows <- incoming[!incoming_keys %in% existing_keys]
  if (nrow(new_rows) && !length(retry_rows) && !schema_migrated) {
    data.table::fwrite(
      new_rows[, names(existing), with = FALSE],
      path,
      sep = "\t",
      quote = FALSE,
      na = "NA",
      append = TRUE,
      col.names = FALSE
    )
  }
  combined <- data.table::rbindlist(
    list(existing, new_rows), use.names = TRUE, fill = FALSE
  )
  data.table::setorder(
    combined, replicate, threshold, method, selection_mode, population
  )
  if (length(retry_rows) || schema_migrated) {
    data.table::fwrite(
      combined, path, sep = "\t", quote = FALSE, na = "NA"
    )
  }
  invisible(combined)
}

merge_hit_regions <- function(hits, radius = 250000L) {
  hits <- unique(data.table::as.data.table(hits)[, .(CHR, BP)])
  hits[, `:=`(
    CHR = as.character(CHR),
    start = pmax(1L, as.integer(BP) - radius),
    end = as.integer(BP) + radius
  )]
  data.table::setorder(hits, CHR, start, end)
  regions <- list()
  index <- 0L
  for (chromosome in unique(hits$CHR)) {
    current <- hits[CHR == chromosome]
    if (!nrow(current)) next
    region_start <- current$start[1L]
    region_end <- current$end[1L]
    for (i in seq_len(nrow(current))[-1L]) {
      if (current$start[i] <= region_end) {
        region_end <- max(region_end, current$end[i])
      } else {
        index <- index + 1L
        regions[[index]] <- data.table::data.table(
          region = paste0("region_", index), CHR = chromosome,
          start = region_start, end = region_end
        )
        region_start <- current$start[i]
        region_end <- current$end[i]
      }
    }
    index <- index + 1L
    regions[[index]] <- data.table::data.table(
      region = paste0("region_", index), CHR = chromosome,
      start = region_start, end = region_end
    )
  }
  data.table::rbindlist(regions)
}
