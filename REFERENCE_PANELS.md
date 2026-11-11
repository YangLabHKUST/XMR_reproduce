# LD reference panels

The real-data workflows use population-specific PLINK `.bed`, `.bim` and `.fam`
files for LD clumping. These genotype panels are separate from the LD-score tables
used for background estimation.

Public GRCh38 1000 Genomes files are available from the
[PLINK 2 resources page](https://www.cog-genomics.org/plink/2.0/resources#1kg_phase3).
[prepare_1kg_reference.sh](prepare_1kg_reference.sh) downloads the genotype,
rsID-annotated variant and corrected sample files, retains autosomal biallelic
variants, splits samples into EUR, EAS, AFR and SAS groups, and applies
`MAF >= 0.005` within each group before writing PLINK binary files.

Install [PLINK 2](https://www.cog-genomics.org/plink/2.0/) and `curl`, then run from
the repository root with an empty output directory:

```bash
bash prepare_1kg_reference.sh negative_controls/1kg_pops
```

The rsID variant download is saved as `all_hg38.pvar.zst`, and the sample download
as `all_hg38.psam`, matching `--pfile all_hg38 vzs`. Keep the PLINK logs and the
downloaded sample metadata. Use [PLINK 1.9](https://www.cog-genomics.org/plink/1.9/)
at `negative_controls/plink_Linux` for the analysis runners.

Expected output prefixes under `negative_controls/1kg_pops/`:

| Population | Prefix | Samples in the study panel | Variants in the study panel |
| --- | --- | ---: | ---: |
| EUR | `EUR_phase3_hg38_autosomes` | 633 | 11,417,608 |
| EAS | `EAS_phase3_hg38_autosomes` | 585 | 10,097,576 |
| AFR | `AFR_phase3_hg38_autosomes` | 893 | 20,914,656 |
| SAS | `SAS_phase3_hg38_autosomes` | 601 | 12,021,508 |

CSA analyses use SAS; BBJ and TPMI use EAS. EUR panels are used for EUR-based
clumping. The sample counts above include related individuals; PLINK 1.9 uses
founders for LD calculations. Preserve pedigree information when preparing panels.

The public downloads may contain updated rsID annotations or pedigree corrections.
This preparation script has not been validated to reproduce the exact historical
panel files or identical clumping results; matching file prefixes or sample counts
alone does not establish equivalence. The four study panels are not included in
the LD-score download.

Simulation LD matrices are separate GRCh37 UK Biobank resources: follow the
[ViLMA instructions](https://github.com/jeffspence/vilma#ld-matrices) and the
[simulation workflow](simulations/README.md).
