#!/bin/bash
# ============================================================================
# Prepare 1000 Genomes Phase 3 (hg38) Reference Panel
# ============================================================================
# This script downloads and processes 1000 Genomes Phase 3 data in GRCh38/hg38
# coordinates, producing per-population PLINK bed/bim/fam files for use as LD
# reference panels.
#
# Prerequisites:
#   - PLINK2 (https://www.cog-genomics.org/plink/2.0/)
#   - wget
#
# Reference:
#   https://cran.r-project.org/web/packages/snpsettest/vignettes/reference_1000Genomes.html
# ============================================================================

# Set your working directory
OUTDIR="./1kg_pops"
mkdir -p "$OUTDIR"
cd "$OUTDIR"

# ============================================================================
# Step 1: Download raw 1000 Genomes files in PLINK2 format
# ============================================================================
# These files are hosted by the PLINK2 developers and contain all 1000 Genomes
# Phase 3 samples genotyped on hg38 coordinates.
#   - all_hg38.pgen.zst : genotype data (zstd-compressed PLINK2 binary)
#   - all_hg38_rs.pvar.zst : variant information (rsIDs, positions, alleles)
#   - hg38_corrected.psam : sample information (sample IDs, sex, population)

wget https://www.dropbox.com/s/j72j6uciq5zuzii/all_hg38.pgen.zst?dl=1 -O all_hg38.pgen.zst
wget https://www.dropbox.com/scl/fi/fn0bcm5oseyuawxfvkcpb/all_hg38_rs.pvar.zst?rlkey=przncwb78rhz4g4ukovocdxaz&dl=1 -O all_hg38_rs.pvar.zst
wget https://www.dropbox.com/s/2e87z6nc4qexjjm/hg38_corrected.psam?dl=1 -O all_hg38.psam

# Also download the sample info files from 1000 Genomes EBI FTP (optional,
# for additional sample metadata)
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_sample_info.xlsx

# ============================================================================
# Step 2: Decompress the genotype file
# ============================================================================
# The .pgen.zst file is compressed with zstandard. PLINK2 has a built-in
# decompressor. This produces the uncompressed .pgen file.

plink2 --zst-decompress all_hg38.pgen.zst > all_hg38.pgen

# ============================================================================
# Step 3: Filter variants and create a clean autosomal dataset
# ============================================================================
# This step:
#   --pfile all_hg38 vzs   : read the pgen/pvar.zst/psam fileset
#   --chr 1-22              : keep only autosomes (chromosomes 1-22)
#   --output-chr 26         : use "chr1" style chromosome naming (UCSC format)
#   --max-alleles 2         : keep only biallelic variants
#   --rm-dup exclude-mismatch : remove duplicate variants with mismatched alleles
#   --set-missing-var-ids   : assign IDs to variants that lack rsIDs,
#                             using format: chr_position_ref_alt
#   --allow-extra-chr       : don't error on non-standard chromosome codes
#   --make-pgen             : output in PLINK2 pgen format

plink2 --pfile all_hg38 vzs \
  --chr 1-22 \
  --output-chr 26 \
  --max-alleles 2 \
  --rm-dup exclude-mismatch \
  --set-missing-var-ids '@_#_$1_$2' \
  --allow-extra-chr \
  --make-pgen \
  --out all_phase3_hg38_autosomes

# ============================================================================
# Step 4: Create per-population sample lists
# ============================================================================
# The .psam file has a super-population column (column 5). We extract sample
# IDs for each of the four major ancestry groups:
#   EUR = European, EAS = East Asian, SAS = South Asian, AFR = African

awk 'NR == 1 || $5 == "EUR" {print $1}' all_hg38.psam > EUR_1kg_samples.txt
awk 'NR == 1 || $5 == "EAS" {print $1}' all_hg38.psam > EAS_1kg_samples.txt
awk 'NR == 1 || $5 == "SAS" {print $1}' all_hg38.psam > SAS_1kg_samples.txt
awk 'NR == 1 || $5 == "AFR" {print $1}' all_hg38.psam > AFR_1kg_samples.txt

# ============================================================================
# Step 5: Subset genotype data by population
# ============================================================================
# For each population, extract its samples from the cleaned autosomal dataset.
# This produces per-population PLINK2 pgen/pvar/psam filesets.

for POP in EUR EAS SAS AFR; do
  plink2 --pfile all_phase3_hg38_autosomes \
    --keep ${POP}_1kg_samples.txt \
    --make-pgen \
    --out ${POP}_phase3_hg38_autosomes
done

# ============================================================================
# Step 6: MAF filter and convert to PLINK1 bed/bim/fam format
# ============================================================================
# For each population, apply a minor allele frequency (MAF) filter of 0.5%
# to remove very rare variants, and convert from PLINK2 pgen format to the
# more widely compatible PLINK1 bed/bim/fam format.
#
# The resulting files (e.g., EUR_phase3_hg38_autosomes.bed/.bim/.fam) are
# used as LD reference panels in downstream analyses.

for POP in EUR EAS SAS AFR; do
  plink2 --pfile ${POP}_phase3_hg38_autosomes \
    --maf 0.005 \
    --make-bed \
    --out ${POP}_phase3_hg38_autosomes
done

echo "Done. Per-population PLINK files are in: $OUTDIR"