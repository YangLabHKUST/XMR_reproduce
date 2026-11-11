#!/usr/bin/env python3
"""Batch greedy LD clumping directly from ViLMA block LD matrices.

The input is a TSV/TSV.GZ table with one row per candidate variant and these
columns: job_id, population, SNP, priority. Smaller priority values are chosen
first. ``population`` may be EUR, EAS, or a comma-separated sequence such as
EUR,EAS; sequential populations require the retained variants to survive
clumping in every listed LD panel.

ViLMA schemas define independent LD blocks. Clumping is therefore performed
within each block, using raw r^2 and the requested physical-distance window.
The SNP reference may be a Step 1/2 table (ID, CHROM) or a Step 5 summary file
(SNP, CHR); both formats are normalized before alignment.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
PARENT_DIR = SCRIPT_DIR.parent
if str(PARENT_DIR) not in sys.path:
    sys.path.insert(0, str(PARENT_DIR))

from step2_compute_vilma_ldscores import (  # noqa: E402
    PreparedSchema,
    load_aligned_block,
    prepare_schema,
)


ROOT = Path("simulations/data")
DEFAULT_EUR_SCHEMA = ROOT / "vilma_ld/wb_hapmap/ld_manifest.txt"
DEFAULT_EAS_SCHEMA = ROOT / "vilma_ld/eas_hapmap/ld_manifest.txt"
REQUIRED_COLUMNS = ["job_id", "population", "SNP", "priority"]
POPULATION_ALIASES = {
    "EUR": "EUR",
    "WB": "EUR",
    "WHITE_BRITISH": "EUR",
    "EAS": "EAS",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Batch LD clumping using ViLMA EUR/EAS block matrices."
    )
    parser.add_argument("--jobs", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--common-snps", type=Path, required=True)
    parser.add_argument("--eur-schema", type=Path, default=DEFAULT_EUR_SCHEMA)
    parser.add_argument("--eas-schema", type=Path, default=DEFAULT_EAS_SCHEMA)
    parser.add_argument("--clump-kb", type=int, default=10_000)
    parser.add_argument("--clump-r2", type=float, default=0.001)
    parser.add_argument("--expected-snps", type=int, default=0)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.clump_kb <= 0:
        parser.error("--clump-kb must be positive")
    if not 0 < args.clump_r2 <= 1:
        parser.error("--clump-r2 must be in (0, 1]")
    if args.expected_snps < 0:
        parser.error("--expected-snps cannot be negative")
    return args


def load_clump_reference(path: Path, expected_snps: int = 0) -> pd.DataFrame:
    """Load and normalize the SNP coordinates used to align ViLMA blocks."""
    if not path.is_file():
        raise FileNotFoundError(f"Clump SNP reference does not exist: {path}")
    reference = pd.read_csv(
        path, sep="\t", compression="infer", low_memory=False
    )
    rename: dict[str, str] = {}
    if "ID" not in reference.columns and "SNP" in reference.columns:
        rename["SNP"] = "ID"
    if "CHROM" not in reference.columns and "CHR" in reference.columns:
        rename["CHR"] = "CHROM"
    reference.rename(columns=rename, inplace=True)

    required = ["ID", "CHROM", "BP", "A1", "A2"]
    missing = [column for column in required if column not in reference.columns]
    if missing:
        raise ValueError(f"Clump SNP reference is missing columns: {missing}")
    optional = [
        column for column in ("EUR_BLOCK", "EAS_BLOCK")
        if column in reference.columns
    ]
    reference = reference[required + optional].copy()
    if expected_snps and len(reference) != expected_snps:
        raise ValueError(
            f"Expected {expected_snps:,} SNPs but found {len(reference):,} in {path}"
        )
    if reference["ID"].isna().any():
        raise ValueError("Clump SNP reference contains missing IDs")
    reference["ID"] = reference["ID"].astype(str)
    duplicated = reference["ID"].duplicated(keep=False)
    if duplicated.any():
        examples = ", ".join(reference.loc[duplicated, "ID"].head(5))
        raise ValueError(
            f"Clump SNP reference contains {duplicated.sum():,} duplicated IDs; "
            f"examples: {examples}"
        )
    positions = pd.to_numeric(reference["BP"], errors="coerce")
    if positions.isna().any() or (positions <= 0).any():
        raise ValueError("Clump SNP reference contains invalid BP values")
    reference["BP"] = positions.astype(np.int64)
    reference["CHROM"] = reference["CHROM"].astype(str)
    reference["A1"] = reference["A1"].astype(str).str.upper()
    reference["A2"] = reference["A2"].astype(str).str.upper()
    return reference


def population_sequence(raw: object) -> tuple[str, ...]:
    text = str(raw).strip().upper().replace(";", ",").replace("+", ",")
    tokens = [token.strip() for token in text.split(",") if token.strip()]
    if not tokens:
        raise ValueError(f"Empty population sequence: {raw!r}")
    result: list[str] = []
    for token in tokens:
        canonical = POPULATION_ALIASES.get(token)
        if canonical is None:
            raise ValueError(
                f"Unsupported population {token!r}; expected EUR/WB or EAS"
            )
        if canonical not in result:
            result.append(canonical)
    return tuple(result)


def validate_jobs(path: Path, common: pd.DataFrame) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Clump jobs file does not exist: {path}")
    jobs = pd.read_csv(path, sep="\t", compression="infer", low_memory=False)
    missing = [column for column in REQUIRED_COLUMNS if column not in jobs.columns]
    if missing:
        raise ValueError(f"Clump jobs file is missing columns: {missing}")
    if jobs.empty:
        raise ValueError("Clump jobs file contains no candidate rows")

    jobs = jobs[REQUIRED_COLUMNS].copy()
    jobs["job_id"] = jobs["job_id"].astype(str)
    jobs["SNP"] = jobs["SNP"].astype(str)
    jobs["population"] = jobs["population"].astype(str)
    jobs["priority"] = pd.to_numeric(jobs["priority"], errors="coerce")
    if (jobs["job_id"].str.len() == 0).any():
        raise ValueError("job_id cannot be empty")
    if not np.isfinite(jobs["priority"].to_numpy(dtype=np.float64)).all():
        raise ValueError("All clump priorities must be finite")

    sequence_by_job: dict[str, str] = {}
    for job_id, group in jobs.groupby("job_id", sort=False):
        sequences = {
            ",".join(population_sequence(value))
            for value in group["population"].unique()
        }
        if len(sequences) != 1:
            raise ValueError(
                f"Job {job_id!r} contains multiple population sequences: "
                f"{sorted(sequences)}"
            )
        sequence_by_job[job_id] = next(iter(sequences))
    jobs["population"] = jobs["job_id"].map(sequence_by_job)

    jobs.sort_values(
        ["job_id", "priority", "SNP"], kind="mergesort", inplace=True
    )
    jobs.drop_duplicates(["job_id", "SNP"], keep="first", inplace=True)

    id_to_index = {snp: index for index, snp in enumerate(common["ID"].astype(str))}
    jobs["global_index"] = jobs["SNP"].map(id_to_index)
    missing_snps = jobs["global_index"].isna()
    if missing_snps.any():
        examples = ", ".join(jobs.loc[missing_snps, "SNP"].head(5))
        raise ValueError(
            f"{missing_snps.sum():,} candidate rows are absent from the common "
            f"SNP list; examples: {examples}"
        )
    jobs["global_index"] = jobs["global_index"].astype(np.int64)
    return jobs.reset_index(drop=True)


def clump_block_candidates(
    candidates: pd.DataFrame,
    block_matrix: np.ndarray,
    local_by_global: dict[int, int],
    bp: np.ndarray,
    window_bp: int,
    r2_threshold: float,
) -> np.ndarray:
    """Return candidate DataFrame indices retained from one job/block."""
    ordered = candidates.sort_values(
        ["priority", "SNP"], kind="mergesort"
    )
    row_indices = ordered.index.to_numpy(dtype=np.int64)
    global_indices = ordered["global_index"].to_numpy(dtype=np.int64)
    local_indices = np.fromiter(
        (local_by_global[int(index)] for index in global_indices),
        dtype=np.int64,
        count=len(global_indices),
    )
    positions = bp[global_indices]
    excluded = np.zeros(len(ordered), dtype=bool)
    retained: list[int] = []

    for lead in range(len(ordered)):
        if excluded[lead]:
            continue
        retained.append(int(row_indices[lead]))
        in_window = np.abs(positions - positions[lead]) <= window_bp
        correlations = block_matrix[local_indices[lead], local_indices]
        linked = np.square(correlations) > r2_threshold
        excluded |= in_window & linked
        excluded[lead] = True
    return np.asarray(retained, dtype=np.int64)


def clump_stage(
    candidates: pd.DataFrame,
    population: str,
    prepared: PreparedSchema,
    common: pd.DataFrame,
    window_bp: int,
    r2_threshold: float,
) -> pd.DataFrame:
    if candidates.empty:
        return candidates.copy()
    working = candidates.copy()
    globals_ = working["global_index"].to_numpy(dtype=np.int64)
    working["block_index"] = prepared.membership[globals_]
    bp = pd.to_numeric(common["BP"], errors="raise").to_numpy(dtype=np.int64)
    retained_indices: list[np.ndarray] = []

    for block_index, block_candidates in working.groupby(
        "block_index", sort=True
    ):
        source_block = prepared.blocks[int(block_index)]
        block_matrix = load_aligned_block(source_block)
        for _, job_candidates in block_candidates.groupby("job_id", sort=False):
            retained_indices.append(
                clump_block_candidates(
                    job_candidates,
                    block_matrix,
                    source_block.local_by_global,
                    bp,
                    window_bp,
                    r2_threshold,
                )
            )
    if not retained_indices:
        return working.iloc[0:0].drop(columns="block_index")
    selected = working.loc[np.concatenate(retained_indices)].copy()
    selected.drop(columns="block_index", inplace=True)
    selected.sort_values(
        ["job_id", "priority", "SNP"], kind="mergesort", inplace=True
    )
    return selected


def run_clumping(
    jobs: pd.DataFrame,
    schemas: dict[str, PreparedSchema],
    common: pd.DataFrame,
    window_bp: int,
    r2_threshold: float,
) -> pd.DataFrame:
    outputs: list[pd.DataFrame] = []
    for sequence_text, sequence_jobs in jobs.groupby("population", sort=False):
        selected = sequence_jobs.copy()
        sequence = tuple(sequence_text.split(","))
        for population in sequence:
            selected = clump_stage(
                selected,
                population,
                schemas[population],
                common,
                window_bp,
                r2_threshold,
            )
        outputs.append(selected)

    selected = pd.concat(outputs, ignore_index=True)
    selected.sort_values(
        ["job_id", "priority", "SNP"], kind="mergesort", inplace=True
    )
    selected["selection_rank"] = (
        selected.groupby("job_id", sort=False).cumcount() + 1
    )
    return selected[
        ["job_id", "population", "SNP", "priority", "selection_rank"]
    ].reset_index(drop=True)


def write_outputs(
    selected: pd.DataFrame,
    jobs: pd.DataFrame,
    output_path: Path,
    summary_path: Path,
    args: argparse.Namespace,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    selected.to_csv(
        output_path,
        sep="\t",
        index=False,
        compression="infer",
        float_format="%.17g",
    )

    candidate_counts = jobs.groupby(
        ["job_id", "population"], sort=False
    ).size().rename("n_candidates")
    selected_counts = selected.groupby(
        ["job_id", "population"], sort=False
    ).size().rename("n_selected")
    summary = candidate_counts.to_frame().join(selected_counts, how="left")
    summary["n_selected"] = summary["n_selected"].fillna(0).astype(int)
    summary.reset_index(inplace=True)
    summary["clump_kb"] = args.clump_kb
    summary["clump_r2"] = args.clump_r2
    summary.to_csv(summary_path, sep="\t", index=False)

    metadata = {
        "jobs": str(args.jobs.resolve()),
        "output": str(output_path.resolve()),
        "summary": str(summary_path.resolve()),
        "common_snps": str(args.common_snps.resolve()),
        "eur_schema": str(args.eur_schema.resolve()),
        "eas_schema": str(args.eas_schema.resolve()),
        "clump_kb": args.clump_kb,
        "clump_r2": args.clump_r2,
        "n_jobs": int(jobs["job_id"].nunique()),
        "n_candidate_rows": int(len(jobs)),
        "n_selected_rows": int(len(selected)),
    }
    metadata_path = summary_path.with_suffix(summary_path.suffix + ".json")
    metadata_path.write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def main() -> int:
    args = parse_args()
    output_path = args.output.expanduser().resolve()
    summary_path = args.summary.expanduser().resolve()
    metadata_path = summary_path.with_suffix(summary_path.suffix + ".json")
    collisions = [
        path for path in (output_path, summary_path, metadata_path) if path.exists()
    ]
    if collisions and not args.overwrite:
        raise FileExistsError(
            "Refusing to overwrite existing clump output(s): "
            + ", ".join(str(path) for path in collisions)
        )

    common = load_clump_reference(
        args.common_snps.expanduser().resolve(), args.expected_snps
    )
    jobs = validate_jobs(args.jobs.expanduser().resolve(), common)
    populations_needed = {
        population
        for sequence in jobs["population"].unique()
        for population in sequence.split(",")
    }
    schemas: dict[str, PreparedSchema] = {}
    if "EUR" in populations_needed:
        schemas["EUR"] = prepare_schema(args.eur_schema, common, "EUR")
    if "EAS" in populations_needed:
        schemas["EAS"] = prepare_schema(args.eas_schema, common, "EAS")

    selected = run_clumping(
        jobs,
        schemas,
        common,
        window_bp=args.clump_kb * 1_000,
        r2_threshold=args.clump_r2,
    )
    write_outputs(selected, jobs, output_path, summary_path, args)
    print(
        f"Clumped {jobs['job_id'].nunique():,} jobs: "
        f"{len(jobs):,} candidates -> {len(selected):,} selected"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
