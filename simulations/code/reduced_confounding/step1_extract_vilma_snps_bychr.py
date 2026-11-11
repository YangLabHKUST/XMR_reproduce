#!/usr/bin/env python3
"""Extract chromosome-specific SNP metadata from ViLMA LD schemas.

The script reads only the schema manifests and their small .var metadata files;
it does not load the LD .npy matrices. By default it extracts chromosomes 2
and 19 from the EAS and White British (named EUR by ViLMA) HapMap3 panels.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Iterable

import pandas as pd


DEFAULT_LD_DIR = Path(
    "simulations/data/vilma_ld"
)
VAR_COLUMNS = ["ID", "CHROM", "BP", "CM", "A1", "A2"]
ALLELE_COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
POPULATION_ALIASES = {
    "eas": ("eas", "east_asian", "east-asian"),
    "eur": (
        "eur",
        "wb",
        "ukb",
        "ukbb",
        "brit",
        "whitebritish",
        "white_british",
        "white-british",
    ),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Extract EAS and EUR/White British variants from ViLMA LD schemas "
            "and construct an allele-aligned common-SNP list."
        )
    )
    parser.add_argument("--ld-dir", type=Path, default=DEFAULT_LD_DIR)
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help=(
            "Output directory. For multiple chromosomes, chromosome-specific "
            "subdirectories are created below it. Default: <ld-dir>/chrN_snps"
        ),
    )
    chromosome_group = parser.add_mutually_exclusive_group()
    chromosome_group.add_argument(
        "--chromosomes",
        nargs="+",
        default=None,
        help=(
            "Chromosomes to extract, separated by spaces or commas "
            "(default: 2 19)"
        ),
    )
    chromosome_group.add_argument(
        "--chromosome",
        action="append",
        dest="single_chromosomes",
        default=None,
        help=(
            "Backward-compatible single-chromosome option; it may be repeated"
        ),
    )
    parser.add_argument(
        "--eas-schema",
        type=Path,
        default=None,
        help="Optional explicit path to the EAS ld_manifest.txt",
    )
    parser.add_argument(
        "--eur-schema",
        type=Path,
        default=None,
        help="Optional explicit path to the EUR/White British ld_manifest.txt",
    )
    parser.add_argument(
        "--drop-ambiguous",
        action="store_true",
        help="Exclude A/T and C/G SNPs from the common-SNP output",
    )
    args = parser.parse_args()
    if args.chromosomes is not None:
        raw_chromosomes = args.chromosomes
    elif args.single_chromosomes is not None:
        raw_chromosomes = args.single_chromosomes
    else:
        raw_chromosomes = ["2", "19"]
    try:
        args.chromosomes = normalized_chromosomes(raw_chromosomes)
    except ValueError as exc:
        parser.error(str(exc))
    del args.single_chromosomes
    return args


def normalized_chromosome(value: object) -> str:
    chrom = str(value).strip().upper()
    chrom = re.sub(r"^CHR", "", chrom)
    if re.fullmatch(r"0*[0-9]+(?:\.0+)?", chrom):
        return str(int(float(chrom)))
    return chrom


def normalized_chromosomes(values: Iterable[object]) -> list[str]:
    chromosomes: list[str] = []
    for value in values:
        for field in str(value).split(","):
            field = field.strip()
            if not field:
                continue
            chromosome = normalized_chromosome(field)
            if not re.fullmatch(r"[A-Z0-9_]+", chromosome):
                raise ValueError(f"Invalid chromosome label: {field!r}")
            if chromosome not in chromosomes:
                chromosomes.append(chromosome)
    if not chromosomes:
        raise ValueError("At least one chromosome must be supplied")
    return chromosomes


def relative_or_absolute(path: Path, base: Path) -> str:
    try:
        return str(path.resolve().relative_to(base.resolve()))
    except ValueError:
        return str(path.resolve())


def schema_entries(schema_path: Path) -> Iterable[tuple[int, Path, Path]]:
    with schema_path.open("r", encoding="utf-8") as handle:
        block_index = 0
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) != 2:
                raise ValueError(
                    f"{schema_path}:{line_number}: expected two fields "
                    "(<variants.var> <ld.npy>)"
                )
            var_path = (schema_path.parent / fields[0]).resolve()
            ld_path = (schema_path.parent / fields[1]).resolve()
            if not var_path.is_file():
                raise FileNotFoundError(f"Variant file not found: {var_path}")
            if not ld_path.is_file():
                raise FileNotFoundError(f"LD matrix file not found: {ld_path}")
            yield block_index, var_path, ld_path
            block_index += 1


def discover_schema(ld_dir: Path, population: str) -> Path:
    manifests = sorted(ld_dir.rglob("ld_manifest.txt"))
    if not manifests:
        raise FileNotFoundError(f"No ld_manifest.txt found below {ld_dir}")

    aliases = POPULATION_ALIASES.get(population.lower(), (population.lower(),))
    patterns = [
        re.compile(rf"(^|[^a-z]){re.escape(alias)}([^a-z]|$)")
        for alias in aliases
    ]
    matches = [
        path
        for path in manifests
        if any(
            pattern.search(str(path.relative_to(ld_dir)).lower())
            for pattern in patterns
        )
    ]
    if len(matches) == 1:
        return matches[0].resolve()
    if not matches:
        available = "\n  ".join(str(x) for x in manifests)
        raise FileNotFoundError(
            f"Could not identify the {population.upper()} schema below {ld_dir}. "
            f"Use --{population}-schema explicitly. Available manifests:\n  "
            f"{available}"
        )
    choices = "\n  ".join(str(x) for x in matches)
    raise RuntimeError(
        f"Multiple {population.upper()} schemas found. Select one with "
        f"--{population}-schema:\n  {choices}"
    )


def read_chromosomes_variants(
    schema_path: Path, chromosomes: Iterable[str]
) -> dict[str, pd.DataFrame]:
    targets = normalized_chromosomes(chromosomes)
    frames: dict[str, list[pd.DataFrame]] = {target: [] for target in targets}

    for block_index, var_path, ld_path in schema_entries(schema_path):
        variants = pd.read_csv(
            var_path,
            sep=r"\s+",
            header=None,
            names=VAR_COLUMNS,
            dtype={"ID": "string", "CHROM": "string", "A1": "string", "A2": "string"},
        )
        normalized = variants["CHROM"].map(normalized_chromosome)
        for target in targets:
            keep = normalized.eq(target)
            if not keep.any():
                continue

            selected = variants.loc[keep].copy()
            selected["CHROM"] = target
            selected["BP"] = pd.to_numeric(
                selected["BP"], errors="raise"
            ).astype("int64")
            selected["A1"] = selected["A1"].str.upper()
            selected["A2"] = selected["A2"].str.upper()
            selected["BLOCK"] = block_index
            selected["VAR_FILE"] = relative_or_absolute(
                var_path, schema_path.parent
            )
            selected["LD_FILE"] = relative_or_absolute(
                ld_path, schema_path.parent
            )
            frames[target].append(selected)

    missing = [target for target in targets if not frames[target]]
    if missing:
        raise RuntimeError(
            f"No variants for chromosome(s) {', '.join(missing)} were found "
            f"in {schema_path}"
        )

    results: dict[str, pd.DataFrame] = {}
    for target in targets:
        result = pd.concat(frames[target], ignore_index=True)
        duplicated = result["ID"].duplicated(keep=False)
        if duplicated.any():
            examples = ", ".join(
                result.loc[duplicated, "ID"].head(5).astype(str)
            )
            raise ValueError(
                f"Duplicate chromosome {target} variant IDs in {schema_path}: "
                f"{duplicated.sum()} rows; examples: {examples}"
            )
        results[target] = result.sort_values(
            ["BP", "ID"], kind="stable"
        ).reset_index(drop=True)
    return results


def read_chromosome_variants(schema_path: Path, chromosome: str) -> pd.DataFrame:
    """Backward-compatible wrapper for callers requesting one chromosome."""
    target = normalized_chromosome(chromosome)
    return read_chromosomes_variants(schema_path, [target])[target]


def complement(allele: str) -> str:
    return "".join(ALLELE_COMPLEMENT.get(base, "N") for base in allele)


def classify_alleles(row: pd.Series) -> str:
    eur = (row["EUR_A1"], row["EUR_A2"])
    eas = (row["EAS_A1"], row["EAS_A2"])
    if eas == eur:
        return "same"
    if eas == eur[::-1]:
        return "swapped"

    eas_complement = (complement(eas[0]), complement(eas[1]))
    if eas_complement == eur:
        return "complement"
    if eas_complement == eur[::-1]:
        return "complement_swapped"
    return "mismatch"


def is_ambiguous(a1: str, a2: str) -> bool:
    return {a1, a2} in ({"A", "T"}, {"C", "G"})


def build_common_variants(
    eas: pd.DataFrame, eur: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    eas_columns = {
        "CHROM": "EAS_CHROM",
        "BP": "EAS_BP",
        "CM": "EAS_CM",
        "A1": "EAS_A1",
        "A2": "EAS_A2",
        "BLOCK": "EAS_BLOCK",
        "VAR_FILE": "EAS_VAR_FILE",
        "LD_FILE": "EAS_LD_FILE",
    }
    eur_columns = {
        "CHROM": "EUR_CHROM",
        "BP": "EUR_BP",
        "CM": "EUR_CM",
        "A1": "EUR_A1",
        "A2": "EUR_A2",
        "BLOCK": "EUR_BLOCK",
        "VAR_FILE": "EUR_VAR_FILE",
        "LD_FILE": "EUR_LD_FILE",
    }
    merged = eur.rename(columns=eur_columns).merge(
        eas.rename(columns=eas_columns), on="ID", how="inner", validate="one_to_one"
    )
    merged["ALLELE_STATUS"] = merged.apply(classify_alleles, axis=1)
    merged["POSITION_MATCH"] = (
        (merged["EUR_CHROM"] == merged["EAS_CHROM"])
        & (merged["EUR_BP"] == merged["EAS_BP"])
    )
    merged["AMBIGUOUS"] = [
        is_ambiguous(a1, a2)
        for a1, a2 in zip(merged["EUR_A1"], merged["EUR_A2"])
    ]
    merged["EAS_SIGN_TO_EUR"] = merged["ALLELE_STATUS"].map(
        {"same": 1, "swapped": -1, "complement": 1, "complement_swapped": -1}
    )

    valid_status = merged["ALLELE_STATUS"].ne("mismatch")
    valid = merged.loc[valid_status & merged["POSITION_MATCH"]].copy()
    rejected = merged.loc[~(valid_status & merged["POSITION_MATCH"])].copy()

    valid.insert(1, "CHROM", valid["EUR_CHROM"])
    valid.insert(2, "BP", valid["EUR_BP"])
    valid.insert(3, "A1", valid["EUR_A1"])
    valid.insert(4, "A2", valid["EUR_A2"])
    valid = valid.sort_values(["BP", "ID"], kind="stable").reset_index(drop=True)
    rejected = rejected.sort_values(["EUR_BP", "ID"], kind="stable").reset_index(drop=True)
    return valid, rejected


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, sep="\t", index=False)
    print(f"Wrote {len(frame):,} rows: {path}")


def main() -> int:
    args = parse_args()
    ld_dir = args.ld_dir.expanduser().resolve()
    if not ld_dir.is_dir():
        raise SystemExit(f"LD directory does not exist: {ld_dir}")

    chromosomes = args.chromosomes
    output_root = (
        args.out_dir.expanduser().resolve()
        if args.out_dir is not None
        else ld_dir
    )

    eas_schema = (
        args.eas_schema.expanduser().resolve()
        if args.eas_schema is not None
        else discover_schema(ld_dir, "eas")
    )
    eur_schema = (
        args.eur_schema.expanduser().resolve()
        if args.eur_schema is not None
        else discover_schema(ld_dir, "eur")
    )
    print(f"EAS schema: {eas_schema}")
    print(f"EUR/White British schema: {eur_schema}")
    print(f"Chromosomes: {', '.join(chromosomes)}")

    eas_by_chromosome = read_chromosomes_variants(eas_schema, chromosomes)
    eur_by_chromosome = read_chromosomes_variants(eur_schema, chromosomes)
    for chromosome in chromosomes:
        if args.out_dir is not None and len(chromosomes) == 1:
            out_dir = output_root
        else:
            out_dir = output_root / f"chr{chromosome}_snps"
        out_dir.mkdir(parents=True, exist_ok=True)

        eas = eas_by_chromosome[chromosome]
        eur = eur_by_chromosome[chromosome]
        common_before_filter, rejected = build_common_variants(eas, eur)
        common = common_before_filter
        if args.drop_ambiguous:
            common = common.loc[~common["AMBIGUOUS"]].copy()

        prefix = f"chr{chromosome}"
        write_tsv(eas, out_dir / f"eas_{prefix}_snps.tsv")
        write_tsv(eur, out_dir / f"eur_{prefix}_snps.tsv")
        write_tsv(common, out_dir / f"eas_eur_{prefix}_common_snps.tsv")
        write_tsv(rejected, out_dir / f"eas_eur_{prefix}_rejected_snps.tsv")

        summary = {
            "chromosome": chromosome,
            "requested_chromosomes": chromosomes,
            "ld_dir": str(ld_dir),
            "output_directory": str(out_dir),
            "eas_schema": str(eas_schema),
            "eur_schema": str(eur_schema),
            "drop_ambiguous": args.drop_ambiguous,
            "eas_snps": int(len(eas)),
            "eur_snps": int(len(eur)),
            "common_by_id": int(len(common_before_filter) + len(rejected)),
            "common_allele_and_position_aligned": int(len(common_before_filter)),
            "ambiguous_in_aligned_common": int(
                common_before_filter["AMBIGUOUS"].sum()
            ),
            "common_written": int(len(common)),
            "rejected": int(len(rejected)),
            "allele_status_counts": {
                str(key): int(value)
                for key, value in pd.concat(
                    [
                        common_before_filter["ALLELE_STATUS"],
                        rejected["ALLELE_STATUS"],
                    ]
                )
                .value_counts()
                .items()
            },
        }
        summary_path = out_dir / f"eas_eur_{prefix}_summary.json"
        summary_path.write_text(
            json.dumps(summary, indent=2) + "\n", encoding="utf-8"
        )
        print(f"Wrote summary: {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
