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
# CAMeRa analysis: SHBG (30830) -> T2D (250.2)
#
# Populations: EUR, EAS, AFR, CSA
# IV thresholds: 5e-5 and 5e-8
# LD clumping: 10000 kb, r2 = 0.001
#
# Formatted input layout:
#   1. Read each population independently, map coordinates with its matched
#      BIM file, and standardise allele direction.
#   2. In each population, select SNPs by that population's SHBG P value
#      and clump them with the matching population LD panel.
#   3. Pool the population-specific hits into regions.
#   4. Within each region, let CAMeRa use four-population fixed-effect
#      meta-analysis (FEMA) to select one shared instrument.
#   5. Extract T2D associations only after the shared IVs are selected,
#      then harmonise and run CAMeRa's cross-population estimator.
#
# Two final-IV modes are analysed at each threshold:
#   A. FEMA_IV: use the regional FEMA instruments directly.
#   B. FEMA_IV_post_population_LD_clump: use FEMA P values to prioritise
#      the FEMA instruments while clumping sequentially with EUR, EAS,
#      AFR and CSA/SAS LD panels. The resulting shared set is independent
#      at 10000 kb and r2 = 0.001 in every panel.
############################################################

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(CAMeRa)
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
  library(ieugwasr)
})

set.seed(1)

# CAMeRa 0.1.0 drops matrix dimensions when a region contains one SNP.
# Keep the official FEMA calculation and preserve a 1 x K matrix.
fixed_effects_meta_analysis_fast_safe <- function(beta_mat, se_mat) {
  if (is.null(dim(beta_mat))) beta_mat <- matrix(beta_mat, nrow = 1L)
  if (is.null(dim(se_mat))) se_mat <- matrix(se_mat, nrow = 1L)
  weights <- 1 / se_mat^2
  beta <- rowSums(beta_mat * weights, na.rm = TRUE) /
    rowSums(weights, na.rm = TRUE)
  se <- sqrt(1 / rowSums(weights, na.rm = TRUE))
  z <- abs(beta / se)
  tibble::tibble(
    nstudy = apply(beta_mat, 1, function(x) sum(is.finite(x))),
    p = pnorm(z, lower.tail = FALSE),
    z = z
  )
}
assignInNamespace(
  "fixed_effects_meta_analysis_fast",
  fixed_effects_meta_analysis_fast_safe,
  ns = "CAMeRa"
)

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
  file.path(.xmr_root, "SHBG_T2D/results/CAMeRa_results_clump10000kb_V2")
)

exposure_code <- "30830"
outcome_code <- "250.2"
pops <- c("EUR", "EAS", "AFR", "CSA")
thresholds <- c(5e-5, 5e-8)
iv_modes <- c("FEMA_IV", "FEMA_IV_post_population_LD_clump")
clump_kb <- 10000L
clump_r2 <- 0.001
region_radius <- 250000L
detected_cores <- parallel::detectCores(logical = FALSE)
if (is.na(detected_cores)) detected_cores <- 1L
mc_cores <- as.integer(get_arg(
  args, "--mc_cores", as.character(min(length(pops), detected_cores))
))
if (is.na(mc_cores) || mc_cores < 1L) {
  stop("--mc_cores must be a positive integer")
}

bfiles <- c(
  EUR = file.path(reference_dir, "1kg_pops/EUR_phase3_hg38_autosomes"),
  EAS = file.path(reference_dir, "1kg_pops/EAS_phase3_hg38_autosomes"),
  AFR = file.path(reference_dir, "1kg_pops/AFR_phase3_hg38_autosomes"),
  CSA = file.path(reference_dir, "1kg_pops/SAS_phase3_hg38_autosomes")
)
plink_bin <- get_arg(args, "--plink_bin", file.path(reference_dir, "plink_Linux"))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
result_file <- file.path(output_dir, "SHBG_T2D_CAMeRa_all_results.tsv")
if (file.exists(result_file)) file.remove(result_file)

find_gwas_file <- function(pop, trait) {
  stem <- file.path(formatted_dir, paste0(pop, "_", trait))
  files <- c(stem, paste0(stem, ".txt"), paste0(stem, ".tsv"), paste0(stem, ".gz"))
  files <- files[file.exists(files)]
  if (!length(files)) stop("Missing GWAS file: ", stem)
  files[1]
}

# BIM files are large, so each population-specific coordinate map is read once.
bim_cache <- new.env(parent = emptyenv())

read_bim <- function(pop) {
  if (exists(pop, bim_cache, inherits = FALSE)) {
    return(get(pop, bim_cache, inherits = FALSE))
  }

  bim <- fread(paste0(bfiles[[pop]], ".bim"), header = FALSE)
  setnames(bim, names(bim)[1:6], c("chr", "SNP", "cm", "position", "ref1", "ref2"))
  bim <- unique(bim[, .(SNP = as.character(SNP), chr, position)], by = "SNP")
  assign(pop, bim, envir = bim_cache)
  bim
}

read_gwas <- function(pop, trait) {
  filename <- find_gwas_file(pop, trait)
  dat <- fread(filename)
  required <- c("SNP", "A1", "A2", "Z", "N", "P")
  missing <- setdiff(required, names(dat))
  if (length(missing)) stop(filename, " is missing: ", paste(missing, collapse = ", "))

  dat[, `:=`(
    SNP = as.character(SNP),
    A1 = toupper(as.character(A1)),
    A2 = toupper(as.character(A2)),
    Z = as.numeric(Z),
    N = as.numeric(N),
    P = as.numeric(P)
  )]
  # Input files have already completed QC; only construct model fields here.
  dat[, `:=`(beta = Z / sqrt(N), se = 1 / sqrt(N))]

  # Use the matched population BIM and standardise alleles alphabetically,
  # exactly as in the EAS and negative-control CAMeRa scripts.
  dat <- merge(dat, read_bim(pop), by = "SNP")
  flip <- dat$A1 > dat$A2
  dat[flip, `:=`(
    beta = -beta,
    Z = -Z,
    A1 = A2,
    A2 = A1
  )]

  dat[, rsid := paste0(chr, ":", position, "_", A1, "_", A2)]
  ans <- dat[, .(
    SNP, rsid, chr, position,
    ea = A1, nea = A2, beta, se, p = P, Z, N,
    eaf = NA_real_, pop = pop,
    id = paste(trait, pop), trait = trait
  )]
  ans[, `:=`(chr = as.character(chr), position = as.integer(position))]
  setindexv(ans, c("chr", "position"))
  setindexv(ans, "rsid")
  ans
}

camera_exp <- setNames(lapply(pops, read_gwas, trait = exposure_code), pops)
camera_out <- setNames(lapply(pops, read_gwas, trait = outcome_code), pops)

failure_row <- function(threshold, iv_mode, error) {
  post_clump <- identical(iv_mode, "FEMA_IV_post_population_LD_clump")
  data.table(
    pops = NA_character_, beta = NA_real_, se = NA_real_, p = NA_real_,
    method = "CAMeRa", exposure = exposure_code, outcome = outcome_code,
    threshold = threshold, iv_mode = iv_mode, nIV = NA_integer_,
    n_fema_iv_before_post_clump = NA_integer_,
    n_pooled_hits = NA_integer_, n_regions = NA_integer_,
    clump_kb = clump_kb, clump_r2 = clump_r2,
    discovery_clump_panel = "population-matched",
    post_fema_clump = post_clump,
    post_fema_clump_priority = if (post_clump) "FEMA P" else "none",
    post_fema_priority_source = if (post_clump) {
      "CAMeRa instrument_fema_regions"
    } else {
      "none"
    },
    post_fema_clump_panels = if (post_clump) {
      "EUR -> EAS -> AFR -> CSA(SAS)"
    } else {
      "none"
    },
    independent_in_all_four_LD_panels = post_clump,
    shared_regional_FEMA_IV = TRUE,
    status = "failed", error = error
  )
}

run_camera <- function(threshold) {
  label <- threshold_label(threshold)
  threshold_dir <- file.path(output_dir, paste0("threshold_", label))
  dir.create(threshold_dir, recursive = TRUE, showWarnings = FALSE)

  # Step 1: select and clump SHBG hits separately in every population.
  hits_by_pop <- setNames(vector("list", length(pops)), pops)
  for (pop in pops) {
    sig <- camera_exp[[pop]][is.finite(p) & p < threshold]
    if (!nrow(sig)) {
      hits_by_pop[[pop]] <- sig
      next
    }

    if (nrow(sig) == 1L) {
      hits <- sig
    } else {
      clumped <- ieugwasr::ld_clump(
        data.frame(rsid = sig$SNP, pval = sig$p),
        clump_kb = clump_kb, clump_r2 = clump_r2, clump_p = threshold,
        bfile = bfiles[[pop]], plink_bin = plink_bin
      )
      hits <- sig[SNP %in% clumped$rsid]
    }
    hits$discovery_pop <- pop
    hits_by_pop[[pop]] <- hits
  }

  all_hits <- rbindlist(hits_by_pop, fill = TRUE)
  pooled_hits <- unique(all_hits[order(p)], by = "SNP")
  if (!nrow(pooled_hits)) stop("No IV remained after population-matched clumping")
  fwrite(
    all_hits,
    file.path(threshold_dir, "CAMeRa_population_specific_P_clumped_hits.tsv"),
    sep = "\t", quote = FALSE
  )

  # Step 2: merge the pooled hits into non-overlapping regional intervals.
  loci <- reduce(GRanges(
    seqnames = pooled_hits$chr,
    ranges = IRanges(
      start = pmax(1L, pooled_hits$position - region_radius),
      end = pooled_hits$position + region_radius
    )
  ))

  # Step 2a: assign exposure SNPs to all loci in one interval-overlap pass per
  # population. This avoids rescanning every full GWAS once for every region.
  loci_dt <- data.table(
    chr = as.character(seqnames(loci)),
    start = as.integer(start(loci)),
    end = as.integer(end(loci)),
    region_name = paste0("region_", seq_along(loci))
  )
  setkey(loci_dt, chr, start, end)

  exposure_by_region <- setNames(lapply(pops, function(pop) {
    points <- copy(camera_exp[[pop]])
    points[, `:=`(
      region_start = as.integer(position),
      region_end = as.integer(position)
    )]
    setkey(points, chr, region_start, region_end)
    overlaps <- foverlaps(
      points, loci_dt,
      by.x = c("chr", "region_start", "region_end"),
      by.y = c("chr", "start", "end"),
      type = "within", nomatch = 0L
    )
    overlaps[, c("start", "end", "region_start", "region_end") := NULL]
    split(overlaps, by = "region_name", keep.by = FALSE, drop = TRUE)
  }), pops)

  # Step 2b: assemble regions in parallel. Outcome tables have an rsid index,
  # so each lookup retrieves only the required SNPs instead of scanning the
  # complete outcome GWAS for every region.
  region_results <- parallel::mclapply(
    loci_dt$region_name,
    function(region_name) {
      exp_region <- setNames(lapply(pops, function(pop) {
        region <- exposure_by_region[[pop]][[region_name]]
        if (is.null(region)) return(camera_exp[[pop]][0])
        unique(region, by = "rsid")
      }), paste("SHBG", pops))

      common_region <- Reduce(intersect, lapply(exp_region, `[[`, "rsid"))
      if (!length(common_region)) return(NULL)

      exp_region <- lapply(exp_region, function(x) {
        ans <- x[rsid %in% common_region]
        setorder(ans, rsid)
        as.data.frame(ans)
      })
      out_region <- setNames(lapply(pops, function(pop) {
        ans <- camera_out[[pop]][
          data.table(rsid = common_region), on = "rsid", nomatch = 0L
        ]
        ans <- unique(ans, by = "rsid")
        setorder(ans, rsid)
        as.data.frame(ans)
      }), paste("T2D", pops))

      # This lead only initialises the object. Regional FEMA reselects the IV.
      exp_all <- rbindlist(exp_region, fill = TRUE)
      lead_vid <- exp_all$rsid[which.max(abs(exp_all$beta / exp_all$se))]
      list(
        region_name = region_name,
        exposure = exp_region,
        outcome = out_region,
        raw_exposure = exp_all[rsid == lead_vid],
        raw_outcome = rbindlist(out_region, fill = TRUE)[rsid == lead_vid]
      )
    },
    mc.cores = min(mc_cores, nrow(loci_dt)),
    mc.preschedule = TRUE
  )
  region_results <- Filter(Negate(is.null), region_results)

  exposure_regions <- setNames(
    lapply(region_results, `[[`, "exposure"),
    vapply(region_results, `[[`, character(1), "region_name")
  )
  outcome_regions <- setNames(
    lapply(region_results, `[[`, "outcome"),
    names(exposure_regions)
  )
  raw_exposure <- setNames(
    lapply(region_results, `[[`, "raw_exposure"),
    names(exposure_regions)
  )
  raw_outcome <- setNames(
    lapply(region_results, `[[`, "raw_outcome"),
    names(exposure_regions)
  )
  if (!length(exposure_regions)) stop("No common regional data for CAMeRa")

  # Step 3: use the official CAMeRa regional fixed-effect meta-analysis.
  fit <- CAMERA$new()
  fit$import_from_local(
    instrument_raw = rbindlist(raw_exposure, fill = TRUE),
    instrument_outcome = rbindlist(raw_outcome, fill = TRUE),
    instrument_regions = exposure_regions,
    instrument_outcome_regions = outcome_regions,
    exposure_ids = paste("SHBG", pops),
    outcome_ids = paste("T2D", pops),
    pops = pops
  )
  fit$fema_regional_instruments(method = "fema")
  if (is.null(fit$instrument_fema) || !nrow(fit$instrument_fema)) {
    stop("CAMeRa regional FEMA returned no instrument")
  }

  fema_instruments <- as.data.table(fit$instrument_fema)
  fwrite(
    fema_instruments,
    file.path(threshold_dir, "CAMeRa_regional_FEMA_shared_IVs.tsv"),
    sep = "\t", quote = FALSE
  )
  saveRDS(fit, file.path(threshold_dir, "CAMeRa_regional_FEMA_base_object.rds"))

  # Recover the official regional FEMA P and Z values. instrument_fema contains
  # population-specific P values, so it must not supply the final-clump priority.
  region_names <- names(fit$instrument_fema_regions)
  if (is.null(region_names)) {
    region_names <- paste0("region_", seq_along(fit$instrument_fema_regions))
  }
  fema_priority <- rbindlist(lapply(seq_along(fit$instrument_fema_regions), function(i) {
    region_dat <- as.data.table(fit$instrument_fema_regions[[i]])
    region_dat <- region_dat[is.finite(z) & is.finite(p)]
    if (!nrow(region_dat)) return(NULL)
    best <- region_dat[which.max(z)]
    data.table(
      region = region_names[i], rsid = as.character(best$rsid[1]),
      p_fema = as.numeric(best$p[1]), z_fema = as.numeric(best$z[1]),
      nstudy = if ("nstudy" %in% names(best)) {
        as.integer(best$nstudy[1])
      } else {
        NA_integer_
      }
    )
  }), fill = TRUE)

  snp_map <- unique(
    fema_instruments[, .(rsid = as.character(rsid), SNP = as.character(SNP))],
    by = "rsid"
  )
  fema_priority <- merge(fema_priority, snp_map, by = "rsid", all.x = TRUE)
  fema_priority <- fema_priority[rsid %in% unique(fema_instruments$rsid)]
  if (!nrow(fema_priority) || anyNA(fema_priority$SNP)) {
    stop("Could not map regional FEMA statistics to original rsIDs")
  }
  fema_priority[, p_fema := pmax(p_fema, .Machine$double.xmin)]
  setorder(fema_priority, p_fema)
  fema_priority <- unique(fema_priority, by = "SNP")
  if (nrow(fema_priority) != uniqueN(fema_instruments$rsid)) {
    stop("Regional FEMA priority table does not match the selected FEMA IV set")
  }
  fwrite(
    fema_priority,
    file.path(threshold_dir, "CAMeRa_regional_FEMA_priority.tsv"),
    sep = "\t", quote = FALSE
  )

  # Run CAMeRa for one final shared-IV set and write mode-specific outputs.
  estimate_mode <- function(instruments, iv_mode, post_clump_panels = "none") {
    mode_dir <- file.path(threshold_dir, iv_mode)
    dir.create(mode_dir, recursive = TRUE, showWarnings = FALSE)

    # Do not require each IV to have exposure/outcome associations in all four
    # populations. CAMeRa harmonises the available associations by population.
    instruments <- as.data.frame(instruments)
    mode_fit <- fit$clone(deep = TRUE)
    mode_fit$instrument_fema <- instruments
    mode_fit$make_outcome_local(exp = mode_fit$instrument_fema)
    mode_fit$harmonise(exp = mode_fit$instrument_fema)
    harmonised <- as.data.table(mode_fit$harmonised_dat)
    if (!nrow(harmonised)) stop("No harmonised SHBG-T2D association remained")

    required_harmonised <- c("pops", "beta.x", "se.x", "beta.y", "se.y")
    missing_harmonised <- setdiff(required_harmonised, names(harmonised))
    if (length(missing_harmonised)) {
      stop("CAMeRa harmonised data lacks: ", paste(missing_harmonised, collapse = ", "))
    }
    harmonised <- harmonised[
      !is.na(pops) & is.finite(beta.x) & is.finite(se.x) & se.x > 0 &
        is.finite(beta.y) & is.finite(se.y) & se.y > 0
    ]
    mode_fit$harmonised_dat <- as.data.frame(harmonised)

    mr <- as.data.table(mode_fit$cross_estimate())
    if (anyNA(mr$Estimate) || anyNA(mr$`Std. Error`) ||
        any(!is.finite(mr$Estimate)) || any(!is.finite(mr$`Std. Error`))) {
      stop("CAMeRa returned a non-finite population-specific estimate")
    }
    fwrite(
      as.data.table(mode_fit$instrument_fema),
      file.path(mode_dir, "CAMeRa_final_shared_IVs.tsv"),
      sep = "\t", quote = FALSE
    )
    fwrite(
      harmonised,
      file.path(mode_dir, "CAMeRa_harmonised_SHBG_T2D.tsv"),
      sep = "\t", quote = FALSE
    )
    saveRDS(mode_fit, file.path(mode_dir, "CAMeRa_SHBG_T2D_object.rds"))

    rename_old <- intersect(
      c("Estimate", "Std. Error", "Pr(>|t|)"), names(mr)
    )
    setnames(
      mr, old = rename_old,
      new = c("beta", "se", "p")[seq_along(rename_old)]
    )

    post_clump <- identical(iv_mode, "FEMA_IV_post_population_LD_clump")
    mr[, `:=`(
      method = "CAMeRa", exposure = exposure_code, outcome = outcome_code,
      threshold = threshold, iv_mode = iv_mode,
      nIV = uniqueN(instruments$rsid),
      n_fema_iv_before_post_clump = uniqueN(fema_instruments$rsid),
      n_harmonised_associations = nrow(harmonised),
      n_pooled_hits = nrow(pooled_hits), n_regions = length(exposure_regions),
      clump_kb = clump_kb, clump_r2 = clump_r2,
      discovery_clump_panel = "population-matched",
      post_fema_clump = post_clump,
      post_fema_clump_priority = if (post_clump) "FEMA P" else "none",
      post_fema_priority_source = if (post_clump) {
        "CAMeRa instrument_fema_regions"
      } else {
        "none"
      },
      post_fema_clump_panels = post_clump_panels,
      independent_in_all_four_LD_panels = post_clump,
      n_clumped_EUR = nrow(hits_by_pop$EUR),
      n_clumped_EAS = nrow(hits_by_pop$EAS),
      n_clumped_AFR = nrow(hits_by_pop$AFR),
      n_clumped_CSA = nrow(hits_by_pop$CSA),
      shared_regional_FEMA_IV = TRUE,
      status = "ok", error = NA_character_
    )]
    mr
  }

  # Mode A: use the official regional FEMA IVs without an extra LD step.
  result_fema <- tryCatch(
    estimate_mode(fema_instruments, "FEMA_IV"),
    error = function(e) failure_row(
      threshold, "FEMA_IV", conditionMessage(e)
    )
  )

  # Mode B: use FEMA P only as priority and sequentially clump the same IV
  # set against all four population-specific panels. clump_p=1 prevents this
  # extra step from imposing a second significance threshold.
  result_post_clump <- tryCatch({
    retained <- copy(fema_priority)
    clump_trace <- list()
    retained_by_panel <- list()

    for (pop in pops) {
      n_before <- nrow(retained)
      if (n_before > 1L) {
        clumped <- ieugwasr::ld_clump(
          data.frame(rsid = retained$SNP, pval = retained$p_fema),
          clump_kb = clump_kb, clump_r2 = clump_r2, clump_p = 1,
          bfile = bfiles[[pop]], plink_bin = plink_bin
        )
        retained <- retained[SNP %in% clumped$rsid]
        setorder(retained, p_fema)
      }

      clump_trace[[pop]] <- data.table(
        panel_population = pop, bfile = bfiles[[pop]],
        n_before = n_before, n_after = nrow(retained),
        n_removed = n_before - nrow(retained)
      )
      retained_by_panel[[pop]] <- copy(retained)[, panel_population := pop]
    }

    fwrite(
      rbindlist(clump_trace),
      file.path(threshold_dir, "CAMeRa_post_FEMA_clump_trace.tsv"),
      sep = "\t", quote = FALSE
    )
    fwrite(
      rbindlist(retained_by_panel, fill = TRUE),
      file.path(threshold_dir, "CAMeRa_post_FEMA_panelwise_retained_IVs.tsv"),
      sep = "\t", quote = FALSE
    )
    if (!nrow(retained)) stop("No FEMA IV remained after four-panel clumping")

    fwrite(
      retained,
      file.path(threshold_dir, "CAMeRa_post_FEMA_final_shared_IVs.tsv"),
      sep = "\t", quote = FALSE
    )
    post_clumped_instruments <- fema_instruments[SNP %in% retained$SNP]
    estimate_mode(
      post_clumped_instruments,
      "FEMA_IV_post_population_LD_clump",
      "EUR -> EAS -> AFR -> CSA(SAS)"
    )
  }, error = function(e) failure_row(
    threshold, "FEMA_IV_post_population_LD_clump", conditionMessage(e)
  ))

  rbindlist(list(result_fema, result_post_clump), fill = TRUE)
}

results <- rbindlist(lapply(thresholds, function(threshold) {
  message("CAMeRa | threshold ", threshold_label(threshold))
  tryCatch(run_camera(threshold), error = function(e) {
    rbindlist(lapply(iv_modes, function(iv_mode) {
      failure_row(threshold, iv_mode, conditionMessage(e))
    }), fill = TRUE)
  })
}), fill = TRUE)

fwrite(results, result_file, sep = "\t", quote = FALSE, na = "NA")
message("CAMeRa finished: ", result_file)
