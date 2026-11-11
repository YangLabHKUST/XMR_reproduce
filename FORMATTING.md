# Formatting GWAS summary statistics

Download the study-specific GWAS files listed in
[GWAS_data_sources.csv](https://github.com/YangLabHKUST/XMR/blob/main/GWAS_data_sources.csv).
Use `MRAPSS::format_data()` after applying the source-specific variant QC and rsID
mapping. Save the formatted tables at the locations in [DATA.md](DATA.md).

The analysis scripts require `SNP`, `A1`, `A2`, `Z`, `N` and `P`.
`A1` is the effect allele and `A2` is the other allele; effect estimates and allele
frequencies must refer to the same effect allele. `P` must be a probability.

For example, for a BBJ quantitative trait:

```r
raw <- data.table::fread("raw_data/GWASsummary.txt.gz")
traits <- read.csv("BBJ_TPMI/raw_data/BBJ/EUR_BBJ_exposure_description.csv")
trait <- traits[as.character(traits$UKB.code) == "21001", ]

formatted <- MRAPSS::format_data(
  raw,
  snp_col = "SNP", b_col = "BETA", se_col = "SE",
  A1_col = "ALLELE1", A2_col = "ALLELE0",
  p_col = "P_BOLT_LMM_INF", freq_col = "A1FREQ",
  n = trait$BBJ.samples, min_freq = 0.01
)
dir.create("BBJ_TPMI/formatted_data", recursive = TRUE, showWarnings = FALSE)
write.table(formatted, "BBJ_TPMI/formatted_data/BBJ_21001",
            quote = FALSE, row.names = FALSE)
```

Adapt the column mapping to the downloaded release:

| Source | Effect / SE | Effect / other allele | P value | Sample size |
| --- | --- | --- | --- | --- |
| BBJ quantitative traits | `BETA` / `SE` | `ALLELE1` / `ALLELE0` | `P_BOLT_LMM_INF` | Trait-specific `BBJ.samples` or `BBJ_sample_size_N` |
| BBJ binary traits and TPMI | `BETA` / `SE` | `Allele2` / `Allele1` | `p.value` | Study- and trait-specific total N |
| Pan-UKB | `beta_<pop>` / `se_<pop>` | `alt` / `ref` | `neglog10_pval_<pop>` (convert as below) | Population- and trait-specific N |
| FinnGen | `beta` / `sebeta` | `alt` / `ref` | `pval` | `FinnGen_sample_size_N` in the auxiliary manifest |
| Neale Lab UKBB quantitative exposures | `beta` / `se` | `alt` / `ref` from `variants.tsv` | `pval` | `n_col = "n_complete_samples"` |

Use `n` for a fixed trait-level sample size or `n_col` for a per-variant sample-size
column. The included [negative-control list](negative_controls/info.csv),
[auxiliary manifest](auxiliary_sensitivity/data/BBJ_PanUKB_FinnGen_68_traits_data_manifest.csv)
and BBJ [exposure](BBJ_TPMI/raw_data/BBJ/EUR_BBJ_exposure_description.csv) /
[outcome](BBJ_TPMI/raw_data/BBJ/EUR_BBJ_outcome_description.csv) descriptions preserve
the population-specific sample sizes and trait identifiers used by the runners.
The BBJ sample sizes are not TPMI sample sizes; obtain TPMI N from its own release
metadata. Preserve phenotype codes such as `1747_5` and `250.2` as identifiers.

For a source column containing `-log10(P)`, pass `log_pval = TRUE` with that column,
or convert it with `P = 10^(-neglog10_P)` before formatting. Convert only once.
For extreme values, the auxiliary Pan-UKB formatter floors converted probabilities
at `.Machine$double.xmin`. Do not pass `-log10(P)` directly as an ordinary P value.

Map variants to rsIDs before formatting when needed. Pan-UKB uses its variant QC
annotation matched on chromosome, position, reference and alternate alleles; Neale
Lab exposures use `variants.tsv` matched by `variant`. TPMI's coordinate-based
mapping uses the hg38 dbSNP153 common-variant annotation. Match the annotation's
genome build to the GWAS coordinates and check duplicate or ambiguous mappings.
PLINK panels are GRCh38, but rsID matching does not require relabelling GRCh37
summary-statistic coordinates as GRCh38.

The formatters use `min_freq = 0.01`. Use the source-specific allele-frequency
column (for example `A1FREQ`, `AF_Allele2` or `af_alt`) and retain the study's QC
choices. Auxiliary Pan-UKB formatting excludes low-confidence EUR variants and
uses case/control sample-size-weighted allele frequencies for binary traits.
These preparation steps precede background estimation and method fitting.
