threshold_label <- function(x) {
  if (isTRUE(all.equal(x, 5e-5))) "5e-05" else "5e-08"
}

scalar_or_na <- function(x) {
  if (is.null(x) || !length(x)) return(NA_real_)
  value <- suppressWarnings(as.numeric(x[1]))
  if (is.finite(value)) value else NA_real_
}

matrix_value <- function(x, i, j) {
  if (is.null(x) || length(dim(x)) != 2L || nrow(x) < i || ncol(x) < j) {
    return(NA_real_)
  }
  scalar_or_na(x[i, j])
}

complement <- function(x) chartr("ACGT", "TGCA", toupper(x))

read_gwas <- function(pop, trait) {
  filename <- file.path(formatted_dir, paste0(pop, "_", trait))
  if (!file.exists(filename)) stop("Missing GWAS file: ", filename)

  dat <- fread(filename)
  required <- c("SNP", "A1", "A2", "Z", "N", "P")
  missing <- setdiff(required, names(dat))
  if (length(missing)) stop(filename, " is missing: ", paste(missing, collapse = ", "))

  dat[, `:=`(
    SNP = as.character(SNP),
    A1 = toupper(as.character(A1)),
    A2 = toupper(as.character(A2)),
    Z = suppressWarnings(as.numeric(Z)),
    N = suppressWarnings(as.numeric(N)),
    P = suppressWarnings(as.numeric(P))
  )]
  # Input files have already completed QC; only construct model fields here.
  dat[, `:=`(se = 1 / sqrt(N), beta = Z / sqrt(N))]
  as.data.frame(dat)
}

load_ld_scores <- function() {
  main_file <- file.path(reference_dir, "LDscore_eas_brit_afr.txt")
  csa_file <- file.path(reference_dir, "LDscore_SAS_EUR.csv")
  if (!file.exists(main_file)) stop("Missing LD-score file: ", main_file)
  if (!file.exists(csa_file)) stop("Missing LD-score file: ", csa_file)

  main <- fread(main_file)
  required_main <- c(
    "rsid", "allele1", "allele2", "EUR", "EAS", "AFR", "EUR_EAS", "AFR_EUR"
  )
  missing_main <- setdiff(required_main, names(main))
  if (length(missing_main)) {
    stop(main_file, " is missing: ", paste(missing_main, collapse = ", "))
  }
  palindromic <- is.na(main$allele1) | is.na(main$allele2) |
    toupper(main$allele1) == complement(main$allele2)
  main <- main[!palindromic]
  main[, `:=`(EUR_AFR = AFR_EUR, EUR_EUR = EUR)]

  csa <- fread(csa_file)
  required_csa <- c("SNP", "EUR", "CSA", "CSA_EUR")
  missing_csa <- setdiff(required_csa, names(csa))
  if (length(missing_csa)) {
    stop(csa_file, " is missing: ", paste(missing_csa, collapse = ", "))
  }
  csa[, EUR_CSA := CSA_EUR]

  list(main = main, csa = csa)
}

get_ldscore <- function(pop) {
  if (pop == "CSA") {
    ans <- ld_scores$csa[, .(
      SNP = as.character(SNP),
      L2.pop1 = as.numeric(EUR),
      L2.pop2 = as.numeric(CSA),
      L12 = as.numeric(EUR_CSA)
    )]
  } else {
    ans <- ld_scores$main[, .(
      SNP = as.character(rsid),
      L2.pop1 = as.numeric(EUR),
      L2.pop2 = as.numeric(get(pop)),
      L12 = as.numeric(get(paste0("EUR_", pop)))
    )]
  }

  ans <- ans[
    nzchar(SNP) &
      is.finite(L2.pop1) & is.finite(L2.pop2) & is.finite(L12)
  ]
  as.data.frame(ans[!duplicated(SNP)])
}

add_two_population_fema <- function(dat) {
  if (!nrow(dat)) return(dat)

  w_eur <- 1 / dat$se.exp.pop1^2
  w_target <- 1 / dat$se.exp.pop2^2
  w_sum <- w_eur + w_target

  dat$b.exp.fema2 <-
    (w_eur * dat$b.exp.pop1 + w_target * dat$b.exp.pop2) / w_sum
  dat$se.exp.fema2 <- sqrt(1 / w_sum)
  dat$z.exp.fema2 <- dat$b.exp.fema2 / dat$se.exp.fema2
  dat$pval.exp.fema2 <- pmax(
    2 * pnorm(abs(dat$z.exp.fema2), lower.tail = FALSE),
    .Machine$double.xmin
  )
  dat
}

prepare_target <- function(pop) {
  exposure_target <- read_gwas(pop, exposure_code)
  outcome_target <- read_gwas(pop, outcome_code)
  exposure_eur <- read_gwas(source_pop, exposure_code)

  aligned_target <- XMR::align_data(
    dat1 = exposure_target,
    dat2 = outcome_target,
    trait1.name = paste0(pop, "_SHBG"),
    trait2.name = paste0(pop, "_T2D"),
    h2.fix.intercept = FALSE,
    LDSC = FALSE
  )$dat
  colnames(aligned_target) <- c(
    "SNP", "A1", "A2", "b.exp.pop2", "b.out.pop2",
    "se.exp.pop2", "se.out.pop2", "pval.exp.pop2", "pval.out.pop2"
  )

  aligned_exposure <- XMR::align_data(
    dat1 = exposure_target,
    dat2 = exposure_eur,
    trait1.name = paste0(pop, "_SHBG"),
    trait2.name = paste0(source_pop, "_SHBG"),
    h2.fix.intercept = FALSE,
    LDSC = FALSE
  )$dat
  colnames(aligned_exposure) <- c(
    "SNP", "A1", "A2", "b.exp.pop2", "b.exp.pop1",
    "se.exp.pop2", "se.exp.pop1", "pval.exp.pop2", "pval.exp.pop1"
  )

  pairwise_data <- merge(
    aligned_target,
    aligned_exposure[, c(
      "SNP", "A1", "A2", "b.exp.pop1", "se.exp.pop1", "pval.exp.pop1"
    )],
    by = c("SNP", "A1", "A2")
  )
  pairwise_data <- merge(pairwise_data, get_ldscore(pop), by = "SNP")
  if (!nrow(pairwise_data)) stop("No EUR-target pairwise SNPs remained after merging")

  reference_candidates <- subset(
    pairwise_data,
    pval.exp.pop1 < max(thresholds)
  )
  if (!nrow(reference_candidates)) {
    stop("No EUR-exposure candidate below the largest IV threshold")
  }

  list(
    data = reference_candidates,
    n_pairwise_merged = nrow(pairwise_data)
  )
}

empty_result <- function(pop, threshold, status,
                         error = NA_character_) {
  data.table(
    method = "XMR",
    selection_mode = "EUR_SELECTED_FEMA_P_CLUMP",
    clump_priority_mode = "FEMA_P",
    reference_population = source_pop,
    population = pop,
    exposure = exposure_code,
    outcome = outcome_code,
    threshold = threshold,
    beta = NA_real_,
    se = NA_real_,
    p = NA_real_,
    tau.sq = NA_real_,
    sigmax11 = NA_real_,
    sigmax12 = NA_real_,
    sigmax22 = NA_real_,
    nIV = NA_real_,
    nvalid = NA_real_,
    IVsignal.sum = NA_real_,
    pi0 = NA_real_,
    n_pairwise_merged = NA_integer_,
    n_candidates = NA_integer_,
    n_clumped_IV = NA_integer_,
    effect_heterogeneity_qc = FALSE,
    candidate_selection = "EUR exposure P below threshold",
    fema_method = "EUR-target signed inverse-variance fixed effect",
    clump_priority = "two-population FEMA P",
    clump_kb = clump_kb,
    clump_r2 = clump_r2,
    clump_panel = paste0(pop, " target panel"),
    dt = 0,
    modified_threshold = threshold,
    status = status,
    error = error
  )
}
