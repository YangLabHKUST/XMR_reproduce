# Data inputs

Obtain the GWAS summary statistics from
[GWAS_data_sources.csv](https://github.com/YangLabHKUST/XMR/blob/main/GWAS_data_sources.csv),
which lists experiment-specific sources and references to the Supplementary tables.
See [FORMATTING.md](FORMATTING.md) for `MRAPSS::format_data()`, sample sizes,
P-value conversion and rsID mapping.
The R scripts read formatted tables with `SNP`, `A1`, `A2`, `Z`, `N` and `P` columns;
`P` is on the probability scale. Additional metadata and panel files are listed below.
The real-data entry points start from these formatted tables.

## Shared real-data references

Place the following resources under `negative_controls/`:

```text
negative_controls/
  plink_Linux
  LDscore_eas_brit_afr.txt
  LDscore_SAS_EUR.csv
  1kg_pops/
    EUR_phase3_hg38_autosomes.{bed,bim,fam}
    EAS_phase3_hg38_autosomes.{bed,bim,fam}
    AFR_phase3_hg38_autosomes.{bed,bim,fam}
    SAS_phase3_hg38_autosomes.{bed,bim,fam}
```

CSA analyses use the SAS reference; BBJ and TPMI use EAS. FinnGen auxiliary analyses use EUR.
The PLINK panels use GRCh38 coordinates. See [REFERENCE_PANELS.md](REFERENCE_PANELS.md)
for the official public source, population splitting and `MAF >= 0.005` preparation.
Place each population's `.bed`, `.bim` and `.fam` files under `1kg_pops/`, with the
prefixes shown above. The public preparation has not been verified to produce
files identical to the study panels.

`LDscore_eas_brit_afr.txt` has columns `chromosome`, `rsid`, `position`, `allele1`, `allele2`,
`AFR`, `EUR`, `EAS`, `AFR_EUR`, `AFR_EAS`, `EUR_EAS`.
`LDscore_SAS_EUR.csv` is tab-delimited despite its extension and has columns
`CHR`, `SNP`, `BP`, `CSA`, `EUR`, `CSA_EUR`, `EUR_CSA`, `EUR_EUR`, `CSA_CSA`.
These score files are distinct resources from the PLINK clumping panels.
Download the two compressed study LD-score files from
[XMR LD-score resources on Zenodo](https://doi.org/10.5281/zenodo.24137958).
Decompress them into `negative_controls/` with their original names.

```bash
# Run in the directory containing the downloaded files.
sha256sum -c SHA256SUMS
gzip -t LDscore_eas_brit_afr.txt.gz LDscore_SAS_EUR.csv.gz
# From the repository root, after placing the .gz files in negative_controls/:
gzip -dk negative_controls/LDscore_eas_brit_afr.txt.gz
gzip -dk negative_controls/LDscore_SAS_EUR.csv.gz
```

On macOS, use `shasum -a 256 -c SHA256SUMS` in place of `sha256sum`.

## Experiment-specific inputs

Four small metadata files are included: the [37-trait negative-control list](negative_controls/info.csv),
the [68-trait auxiliary manifest](auxiliary_sensitivity/data/BBJ_PanUKB_FinnGen_68_traits_data_manifest.csv),
and the BBJ [7-exposure](BBJ_TPMI/raw_data/BBJ/EUR_BBJ_exposure_description.csv) and
[15-outcome](BBJ_TPMI/raw_data/BBJ/EUR_BBJ_outcome_description.csv) descriptions.
Their identifiers and row order are retained for the analysis scripts. Keep these
filenames unchanged. The auxiliary manifest describes the candidate traits; the
final comparison uses the available, matched results selected by the plotting code.

| Analysis | Required files relative to the repository root |
| --- | --- |
| Negative controls | `negative_controls/info.csv`; `negative_controls/formatted_data/{EUR,AFR,CSA}_<phenocode>` |
| Auxiliary sensitivity | `auxiliary_sensitivity/data/BBJ_PanUKB_FinnGen_68_traits_data_manifest.csv`; `auxiliary_sensitivity/data/formatted/{Fin,PanUKB,BBJ}_<trait>` |
| SHBG–T2D | `SHBG_T2D/formatted_data/{EUR,EAS,AFR,CSA}_{30830,250.2}` |
| BBJ–TPMI | `BBJ_TPMI/formatted_data/{EUR,BBJ,TW}_<phenocode>`; `BBJ_TPMI/raw_data/BBJ/EUR_BBJ_exposure_description.csv`; `BBJ_TPMI/raw_data/BBJ/EUR_BBJ_outcome_description.csv` |
| BBJ–TPMI source-outcome comparisons | Pan-UKB formatted outcomes under `auxiliary_sensitivity/data/formatted/`, as selected in the runner configuration |

## Generated background parameters

Background parameters are intermediate outputs. They do not need to be supplied before
running an experiment from its formatted inputs.

| Analysis | How backgrounds are generated |
| --- | --- |
| Negative controls | `code/prepare_background.R` writes MRBEE/TL target tables to `results/background_MRBEE_TL_target/`. The method runner estimates its EUR source backgrounds. XMR and EUR MR-APSS estimate their own backgrounds within their runners. |
| Auxiliary sensitivity | XMR and MRBEE/TL estimate their required backgrounds within their runners. |
| SHBG–T2D | `code/prepare_background.R` writes comparator `C` and `Omega0` to `results/SHBG-T2D-bgparas.RData`. The XMR pipeline estimates its own zero-cross-intercept backgrounds. |
| BBJ–TPMI | XMR, MRBEE/TL and MR-APSS estimate their backgrounds within their runners. |
| Simulations | The launcher runs the background stage before method fitting and passes the output location to each method. |

Paths in this table are relative to the corresponding analysis folder. The real-data
`code/run_analysis.R` entry points execute background preparation before dependent methods.

## Simulation LD matrices

Download the White British and EAS HapMap3 reference matrices from
[ViLMA LD matrices](https://github.com/jeffspence/vilma#ld-matrices). These are UK Biobank
references in GRCh37 coordinates, separate from the real-data PLINK panels.
Keep each manifest together with its referenced `.npy` and `.var` block files:

```text
simulations/data/vilma_ld/
  wb_hapmap/ld_manifest.txt
  eas_hapmap/ld_manifest.txt
  chr2_snps/eas_eur_chr2_ldscores.tsv.gz
  chr19_snps/eas_eur_chr19_ldscores.tsv.gz
```

The Python preparation steps create the chromosome-specific SNP and score tables from the
reference blocks. Chromosomes 2 and 19 contain 69,151 and 15,423 aligned common SNPs,
respectively (84,574 total). The simulation launcher accepts an alternative data root through
`XMR_SIM_DATA_ROOT`.
