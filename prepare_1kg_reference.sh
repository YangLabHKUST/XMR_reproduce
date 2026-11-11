#!/usr/bin/env bash
set -euo pipefail

# Public 1000 Genomes GRCh38 panels; see REFERENCE_PANELS.md.
# Usage from the repository root: bash prepare_1kg_reference.sh
# Optional first argument: an empty output directory.
outdir=${1:-negative_controls/1kg_pops}
plink2_bin=${PLINK2_BIN:-plink2}
plink2_bin=$(command -v "$plink2_bin")
case "$plink2_bin" in
  /*) ;;
  *) plink2_bin="$(cd "$(dirname "$plink2_bin")" && pwd)/$(basename "$plink2_bin")" ;;
esac
command -v curl >/dev/null
mkdir -p "$outdir"
cd "$outdir"
if [[ -n $(ls -A .) ]]; then
  printf '%s\n' 'Use an empty output directory to avoid replacing an existing panel.' >&2
  exit 1
fi

download() {
  curl --fail --location --retry 3 --output "$2.partial" "$1"
  test -s "$2.partial"
  mv "$2.partial" "$2"
}

# Download URLs are listed on the PLINK 2 resources page.
download 'https://www.dropbox.com/s/j72j6uciq5zuzii/all_hg38.pgen.zst?dl=1' all_hg38.pgen.zst
download 'https://www.dropbox.com/scl/fi/fn0bcm5oseyuawxfvkcpb/all_hg38_rs.pvar.zst?dl=1&rlkey=przncwb78rhz4g4ukovocdxaz' all_hg38.pvar.zst
download 'https://www.dropbox.com/scl/fi/u5udzzaibgyvxzfnjcvjc/hg38_corrected.psam?dl=1&rlkey=oecjnk4vmbhc8b1p202l0ih4x' all_hg38.psam

"$plink2_bin" --zst-decompress all_hg38.pgen.zst > all_hg38.pgen.partial
test -s all_hg38.pgen.partial
mv all_hg38.pgen.partial all_hg38.pgen

"$plink2_bin" --pfile all_hg38 vzs \
  --chr 1-22 --output-chr 26 --max-alleles 2 \
  --rm-dup exclude-mismatch --set-missing-var-ids '@_#_$1_$2' \
  --allow-extra-chr --make-pgen --out all_phase3_hg38_autosomes

for pop in EUR EAS AFR SAS; do
  # Read named columns so changes in sample-column order cannot change membership.
  awk -v pop="$pop" '
    NR == 1 {
      for (i = 1; i <= NF; i++) {
        name = $i; sub(/^#/, "", name)
        if (name == "IID") iid = i
        if (name == "FID") fid = i
        if (name == "SuperPop") superpop = i
      }
      if (!iid || !superpop) {
        print "Expected IID and SuperPop columns in the PSAM file" > "/dev/stderr"
        exit 1
      }
      print fid ? "#FID\tIID" : "#IID"
      next
    }
    $superpop == pop {
      if (fid) print $fid "\t" $iid; else print $iid
      count++
    }
    END { if (!count) exit 1 }
  ' all_hg38.psam > "${pop}_1kg_samples.txt"

  "$plink2_bin" --pfile all_phase3_hg38_autosomes \
    --keep "${pop}_1kg_samples.txt" --make-pgen \
    --out "${pop}_phase3_hg38_autosomes"
  "$plink2_bin" --pfile "${pop}_phase3_hg38_autosomes" \
    --maf 0.005 --make-bed --out "${pop}_phase3_hg38_autosomes"
done

printf '%s\n' 'Created EUR, EAS, AFR and SAS PLINK panels in the output directory.'
