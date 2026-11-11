#!/usr/bin/env python3
"""Create LD-independent SNP lists from ViLMA LD matrices.

The script consumes the chromosome-specific LD-score files produced by
``step2_compute_vilma_ldscores.py``. By default it creates separate EUR-only
and EAS-only independent-SNP lists plus a joint EUR/EAS list. Population-only
lists are ranked by their own LD score; the joint list is ranked by EUR LD
score. Every Step 2 LD-score column is retained in each output.

ViLMA manifests define independent LD blocks. Correlations between different
blocks, and between chromosomes, are therefore treated as zero.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from step2_compute_vilma_ldscores import (
    PreparedSchema,
    derived_summary_path,
    file_sha256,
    load_aligned_block,
    load_common_snps,
    normalized_chromosome,
    normalized_chromosomes,
    prepare_schema,
)


DEFAULT_LD_DIR = Path(
    "simulations/data/vilma_ld"
)
DEFAULT_CLUMP_R2 = 0.001
DEFAULT_CLUMP_KB = 10_000
REQUIRED_LDSCORE_COLUMNS = ["L2_EUR", "L2_EAS", "L12_EUR_EAS"]
POPULATION_ALIASES = {
    "EUR": "EUR",
    "WB": "EUR",
    "WHITE_BRITISH": "EUR",
    "WHITE-BRITISH": "EUR",
    "EAS": "EAS",
}


@dataclass(frozen=True)
class ChromosomeJob:
    chromosome: str
    ldscore_path: Path
    output_path: Path
    summary_path: Path


def parse_populations(values: Iterable[object]) -> list[str]:
    populations: list[str] = []
    for value in values:
        text = str(value).strip().upper().replace(";", ",").replace("+", ",")
        for token in text.split(","):
            token = token.strip()
            if not token:
                continue
            canonical = POPULATION_ALIASES.get(token)
            if canonical is None:
                raise ValueError(
                    f"Unsupported clumping population {token!r}; "
                    "use EUR/WB and/or EAS"
                )
            if canonical not in populations:
                populations.append(canonical)
    if not populations:
        raise ValueError("At least one clumping population must be supplied")
    return populations


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create chromosome-specific and combined LD-independent SNP "
            "lists from ViLMA matrices while retaining Step 2 LD scores."
        )
    )
    parser.add_argument("--ld-dir", type=Path, default=DEFAULT_LD_DIR)
    parser.add_argument(
        "--chromosomes",
        nargs="+",
        default=["2", "19"],
        help="Chromosomes separated by spaces or commas (default: 2 19)",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help=(
            "Optional output root. Per-chromosome chrN_snps directories and "
            "the combined file are written below it. By default, chromosome "
            "files are written beside Step 2 inputs and the combined file "
            "below <ld-dir>/independent_snps."
        ),
    )
    parser.add_argument(
        "--combined-output",
        type=Path,
        default=None,
        help="Optional explicit path for the combined independent-SNP list",
    )
    parser.add_argument(
        "--eur-schema",
        type=Path,
        default=None,
        help="EUR/White British manifest (default: <ld-dir>/wb_hapmap/ld_manifest.txt)",
    )
    parser.add_argument(
        "--eas-schema",
        type=Path,
        default=None,
        help="EAS manifest (default: <ld-dir>/eas_hapmap/ld_manifest.txt)",
    )
    parser.add_argument(
        "--clump-populations",
        nargs="+",
        default=None,
        help=(
            "Run only this panel combination instead of the default EUR-only, "
            "EAS-only, and joint EUR/EAS lists"
        ),
    )
    parser.add_argument(
        "--clump-r2",
        type=float,
        default=DEFAULT_CLUMP_R2,
        help="Require raw pairwise r-squared strictly below this value (default: 0.001)",
    )
    parser.add_argument(
        "--clump-kb",
        type=int,
        default=DEFAULT_CLUMP_KB,
        help="Symmetric physical clumping window in kb (default: 10000 = 10 Mb)",
    )
    parser.add_argument(
        "--priority-column",
        default=None,
        help=(
            "Step 2 column used to rank leads (default: L2_EUR for EUR/joint, "
            "L2_EAS for EAS-only)"
        ),
    )
    parser.add_argument(
        "--priority-order",
        choices=("descending", "ascending"),
        default="descending",
        help="Lead-SNP ordering for the priority column (default: descending)",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=20,
        help="Print progress after this many LD blocks; use 0 to disable",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace existing Step 3 outputs",
    )
    args = parser.parse_args()

    try:
        args.chromosomes = normalized_chromosomes(args.chromosomes)
        if args.clump_populations is not None:
            args.clump_populations = parse_populations(args.clump_populations)
    except ValueError as exc:
        parser.error(str(exc))
    if not 0 < args.clump_r2 <= 1:
        parser.error("--clump-r2 must be in (0, 1]")
    if args.clump_kb <= 0:
        parser.error("--clump-kb must be positive")
    if args.progress_every < 0:
        parser.error("--progress-every cannot be negative")
    return args


def resolve_jobs(
    ld_dir: Path,
    output_root: Path | None,
    chromosomes: list[str],
    populations: list[str],
) -> list[ChromosomeJob]:
    tag = population_tag(populations)
    jobs: list[ChromosomeJob] = []
    for chromosome in chromosomes:
        input_dir = ld_dir / f"chr{chromosome}_snps"
        ldscore_path = (
            input_dir / f"eas_eur_chr{chromosome}_ldscores.tsv.gz"
        ).resolve()
        output_dir = (
            input_dir
            if output_root is None
            else output_root / f"chr{chromosome}_snps"
        )
        output_path = (
            output_dir / f"{tag}_chr{chromosome}_independent_snps.tsv.gz"
        ).resolve()
        jobs.append(
            ChromosomeJob(
                chromosome=chromosome,
                ldscore_path=ldscore_path,
                output_path=output_path,
                summary_path=derived_summary_path(output_path),
            )
        )
    return jobs


def default_combined_output(
    ld_dir: Path,
    output_root: Path | None,
    chromosomes: list[str],
    populations: list[str],
) -> Path:
    chromosome_label = "_".join(f"chr{chromosome}" for chromosome in chromosomes)
    tag = population_tag(populations)
    directory = (
        ld_dir / "independent_snps"
        if output_root is None
        else output_root
    )
    return (directory / f"{tag}_{chromosome_label}_independent_snps.tsv.gz").resolve()


def population_tag(populations: list[str]) -> str:
    if populations == ["EUR"]:
        return "eur"
    if populations == ["EAS"]:
        return "eas"
    if set(populations) == {"EUR", "EAS"}:
        return "eas_eur"
    raise ValueError(f"Unsupported population combination: {populations}")


def preflight(
    jobs: list[ChromosomeJob],
    combined_output: Path,
    overwrite: bool,
) -> bool:
    missing = [job.ldscore_path for job in jobs if not job.ldscore_path.is_file()]
    if missing:
        text = "\n  ".join(str(path) for path in missing)
        raise FileNotFoundError(f"Missing Step 2 LD-score file(s):\n  {text}")

    if overwrite:
        return False
    output_pairs = [
        pair
        for job in jobs
        for pair in [(job.output_path, job.summary_path)]
    ]
    output_pairs.append((combined_output, derived_summary_path(combined_output)))
    states = [path.exists() for pair in output_pairs for path in pair]
    if all(states):
        return True
    if any(states):
        existing = [
            path
            for pair in output_pairs
            for path in pair
            if path.exists()
        ]
        text = "\n  ".join(str(path) for path in existing)
        raise FileExistsError(
            f"Incomplete Step 3 output set:\n  {text}\n"
            "Complete or move the existing files, or use --overwrite."
        )
    return False


def load_ldscores(
    path: Path,
    chromosome: str,
    priority_column: str,
) -> pd.DataFrame:
    scores = load_common_snps(path, expected_snps=0)
    if scores.empty:
        raise ValueError(f"Step 2 LD-score file is empty: {path}")
    missing = [
        column
        for column in REQUIRED_LDSCORE_COLUMNS
        if column not in scores.columns
    ]
    if missing:
        raise ValueError(f"{path} is missing LD-score columns: {missing}")
    if priority_column not in scores.columns:
        raise ValueError(
            f"Priority column {priority_column!r} is absent from {path}"
        )

    observed = {
        normalized_chromosome(value)
        for value in scores["CHROM"].dropna().unique()
    }
    if observed != {chromosome}:
        raise ValueError(
            f"{path} contains chromosome labels {sorted(observed)}; "
            f"expected only chromosome {chromosome}"
        )

    numeric_columns = list(REQUIRED_LDSCORE_COLUMNS)
    numeric_columns.extend(
        column
        for column in (
            "L2_EUR_RAW",
            "L2_EAS_RAW",
            "L12_EUR_EAS_RAW",
            "L2.pop1",
            "L2.pop2",
            "L12",
        )
        if column in scores.columns
    )
    if priority_column not in numeric_columns:
        numeric_columns.append(priority_column)
    for column in numeric_columns:
        values = pd.to_numeric(scores[column], errors="coerce").to_numpy(
            dtype=np.float64
        )
        if not np.isfinite(values).all():
            raise ValueError(f"Non-finite values found in {path}:{column}")
        scores[column] = values

    bp = pd.to_numeric(scores["BP"], errors="coerce").to_numpy(dtype=np.float64)
    if not np.isfinite(bp).all() or not np.equal(bp, np.floor(bp)).all():
        raise ValueError(f"BP must contain finite integer positions in {path}")
    scores["BP"] = bp.astype(np.int64)
    return scores


def priority_ranks(
    scores: pd.DataFrame,
    priority_column: str,
    priority_order: str,
) -> np.ndarray:
    ordering = pd.DataFrame(
        {
            "global_index": np.arange(len(scores), dtype=np.int64),
            "priority": scores[priority_column].to_numpy(dtype=np.float64),
            "BP": scores["BP"].to_numpy(dtype=np.int64),
            "ID": scores["ID"].astype(str).to_numpy(),
        }
    )
    ordering.sort_values(
        ["priority", "BP", "ID"],
        ascending=[priority_order == "ascending", True, True],
        kind="mergesort",
        inplace=True,
    )
    ranks = np.empty(len(scores), dtype=np.int64)
    ranks[ordering["global_index"].to_numpy(dtype=np.int64)] = np.arange(
        len(scores), dtype=np.int64
    )
    return ranks


def clump_block(
    candidates: np.ndarray,
    block_matrix: np.ndarray,
    local_by_global: dict[int, int],
    bp: np.ndarray,
    ranks: np.ndarray,
    window_bp: int,
    r2_threshold: float,
) -> np.ndarray:
    ordered = candidates[
        np.argsort(ranks[candidates], kind="stable")
    ]
    local = np.fromiter(
        (local_by_global[int(index)] for index in ordered),
        dtype=np.int64,
        count=len(ordered),
    )
    positions = bp[ordered]
    excluded = np.zeros(len(ordered), dtype=bool)
    retained: list[int] = []

    for lead in range(len(ordered)):
        if excluded[lead]:
            continue
        retained.append(int(ordered[lead]))
        in_window = np.abs(positions - positions[lead]) <= window_bp
        linked = np.square(block_matrix[local[lead], local]) >= r2_threshold
        excluded |= in_window & linked
        excluded[lead] = True
    return np.asarray(retained, dtype=np.int64)


def clump_stage(
    selected: np.ndarray,
    prepared: PreparedSchema,
    bp: np.ndarray,
    ranks: np.ndarray,
    window_bp: int,
    r2_threshold: float,
    progress_every: int,
) -> np.ndarray:
    memberships = prepared.membership[selected]
    block_indices = np.unique(memberships)
    retained: list[np.ndarray] = []

    for count, block_index in enumerate(block_indices, start=1):
        candidates = selected[memberships == block_index]
        block = prepared.blocks[int(block_index)]
        retained.append(
            clump_block(
                candidates,
                load_aligned_block(block),
                block.local_by_global,
                bp,
                ranks,
                window_bp,
                r2_threshold,
            )
        )
        if progress_every and (
            count % progress_every == 0 or count == len(block_indices)
        ):
            print(
                f"{prepared.population} clumping: processed "
                f"{count:,}/{len(block_indices):,} blocks"
            )

    if not retained:
        return np.empty(0, dtype=np.int64)
    output = np.concatenate(retained)
    return output[np.argsort(ranks[output], kind="stable")]


def verify_independence(
    selected: np.ndarray,
    prepared: PreparedSchema,
    bp: np.ndarray,
    window_bp: int,
    r2_threshold: float,
) -> dict[str, object]:
    memberships = prepared.membership[selected]
    maximum_r2 = 0.0
    pairs_in_window = 0
    violating_pairs = 0

    for block_index in np.unique(memberships):
        candidates = selected[memberships == block_index]
        if len(candidates) < 2:
            continue
        block = prepared.blocks[int(block_index)]
        matrix = load_aligned_block(block)
        local = np.fromiter(
            (block.local_by_global[int(index)] for index in candidates),
            dtype=np.int64,
            count=len(candidates),
        )
        positions = bp[candidates]
        for lead in range(len(candidates) - 1):
            later = np.arange(lead + 1, len(candidates), dtype=np.int64)
            in_window = np.abs(positions[later] - positions[lead]) <= window_bp
            if not in_window.any():
                continue
            values = np.square(matrix[local[lead], local[later[in_window]]])
            pairs_in_window += int(values.size)
            maximum_r2 = max(maximum_r2, float(np.max(values)))
            violating_pairs += int(np.count_nonzero(values >= r2_threshold))

    if violating_pairs:
        raise AssertionError(
            f"{prepared.population} verification found {violating_pairs:,} "
            f"retained pairs with r^2 >= {r2_threshold:g} within "
            f"{window_bp / 1_000_000:g} Mb"
        )
    return {
        "population": prepared.population,
        "pairs_checked_within_same_block_and_window": pairs_in_window,
        "violating_pairs": violating_pairs,
        "maximum_retained_pair_r2_within_same_block_and_window": maximum_r2,
        "cross_block_correlations_assumed_zero": True,
    }


def ordered_output_columns(frame: pd.DataFrame) -> list[str]:
    preferred = [
        "ID",
        "CHROM",
        "BP",
        "A1",
        "A2",
        "AMBIGUOUS",
        "L2_EUR",
        "L2_EAS",
        "L12_EUR_EAS",
        "L2_EUR_RAW",
        "L2_EAS_RAW",
        "L12_EUR_EAS_RAW",
        "L2.pop1",
        "L2.pop2",
        "L12",
        "CLUMP_PRIORITY",
        "CLUMP_RANK",
        "CLUMP_POPULATIONS",
        "CLUMP_R2_THRESHOLD",
        "CLUMP_WINDOW_KB",
    ]
    columns = [column for column in preferred if column in frame.columns]
    columns.extend(column for column in frame.columns if column not in columns)
    return columns


def atomic_write_table(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    frame.to_csv(
        temporary,
        sep="\t",
        index=False,
        float_format="%.12g",
        compression="gzip" if str(path).endswith(".gz") else None,
    )
    temporary.replace(path)


def atomic_write_json(data: dict[str, object], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(data, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def process_chromosome(
    job: ChromosomeJob,
    schema_paths: dict[str, Path],
    args: argparse.Namespace,
) -> tuple[pd.DataFrame, dict[str, object]]:
    scores = load_ldscores(
        job.ldscore_path,
        job.chromosome,
        args.priority_column,
    )
    ranks = priority_ranks(
        scores,
        args.priority_column,
        args.priority_order,
    )
    bp = scores["BP"].to_numpy(dtype=np.int64)
    schemas = {
        population: prepare_schema(schema_paths[population], scores, population)
        for population in args.clump_populations
    }

    selected = np.arange(len(scores), dtype=np.int64)
    stage_counts: list[dict[str, object]] = []
    window_bp = args.clump_kb * 1_000
    print(
        f"Chromosome {job.chromosome}: clumping {len(selected):,} SNPs; "
        f"priority={args.priority_column} {args.priority_order}"
    )
    for population in args.clump_populations:
        before = len(selected)
        selected = clump_stage(
            selected,
            schemas[population],
            bp,
            ranks,
            window_bp,
            args.clump_r2,
            args.progress_every,
        )
        stage_counts.append(
            {
                "population": population,
                "num_input": before,
                "num_retained": int(len(selected)),
            }
        )

    verification = [
        verify_independence(
            selected,
            schemas[population],
            bp,
            window_bp,
            args.clump_r2,
        )
        for population in args.clump_populations
    ]
    selected = selected[np.argsort(ranks[selected], kind="stable")]
    output = scores.iloc[selected].copy()
    output["CLUMP_PRIORITY"] = output[args.priority_column].to_numpy(
        dtype=np.float64
    )
    output["CLUMP_RANK"] = np.arange(1, len(output) + 1, dtype=np.int64)
    output["CLUMP_POPULATIONS"] = ",".join(args.clump_populations)
    output["CLUMP_R2_THRESHOLD"] = args.clump_r2
    output["CLUMP_WINDOW_KB"] = args.clump_kb
    for population in args.clump_populations:
        output[f"{population}_LD_BLOCK"] = schemas[population].membership[selected]
    output = output.loc[:, ordered_output_columns(output)]

    atomic_write_table(output, job.output_path)
    summary: dict[str, object] = {
        "chromosome": job.chromosome,
        "input_ldscore_file": str(job.ldscore_path),
        "input_ldscore_file_sha256": file_sha256(job.ldscore_path),
        "output_file": str(job.output_path),
        "num_candidates": int(len(scores)),
        "num_independent_snps": int(len(output)),
        "retained_fraction": float(len(output) / len(scores)),
        "clump_r2": float(args.clump_r2),
        "clump_window_kb": int(args.clump_kb),
        "clump_window_bp": int(window_bp),
        "clump_populations": args.clump_populations,
        "priority_column": args.priority_column,
        "priority_order": args.priority_order,
        "tie_breakers": ["BP ascending", "ID ascending"],
        "ld_for_clumping": (
            "raw squared correlations reconstructed from ViLMA LD matrices; "
            "finite-reference corrected LD scores are not used as pairwise r^2"
        ),
        "retention_rule": (
            "retain pairs only when raw r^2 is strictly below clump_r2 "
            "within the physical window"
        ),
        "step2_ldscore_columns_retained": [
            column
            for column in scores.columns
            if column.startswith("L2") or column.startswith("L12")
        ],
        "stage_counts": stage_counts,
        "verification": verification,
        "schemas": {
            population: {
                "path": str(schema_paths[population]),
                "sha256": file_sha256(schema_paths[population]),
                "num_blocks_used": int(len(schemas[population].blocks)),
                "orientation_counts": schemas[population].orientation_counts,
            }
            for population in args.clump_populations
        },
    }
    atomic_write_json(summary, job.summary_path)
    print(
        f"Chromosome {job.chromosome}: retained {len(output):,}/"
        f"{len(scores):,} SNPs -> {job.output_path}"
    )
    return output, summary


def write_combined_output(
    outputs: list[pd.DataFrame],
    summaries: list[dict[str, object]],
    combined_output: Path,
    args: argparse.Namespace,
) -> None:
    combined = pd.concat(outputs, ignore_index=True)
    combined.rename(
        columns={"CLUMP_RANK": "CHROMOSOME_CLUMP_RANK"}, inplace=True
    )
    chromosome_order = {
        chromosome: index for index, chromosome in enumerate(args.chromosomes)
    }
    combined["_CHROMOSOME_ORDER"] = [
        chromosome_order[normalized_chromosome(value)]
        for value in combined["CHROM"]
    ]
    combined.sort_values(
        ["CLUMP_PRIORITY", "_CHROMOSOME_ORDER", "BP", "ID"],
        ascending=[args.priority_order == "ascending", True, True, True],
        kind="mergesort",
        inplace=True,
    )
    combined.drop(columns="_CHROMOSOME_ORDER", inplace=True)
    combined.reset_index(drop=True, inplace=True)
    combined["CLUMP_RANK"] = np.arange(1, len(combined) + 1, dtype=np.int64)
    combined["COMBINED_LIST_INDEX"] = combined["CLUMP_RANK"]
    combined = combined.loc[:, ordered_output_columns(combined)]
    atomic_write_table(combined, combined_output)
    combined_summary = {
        "requested_chromosomes": args.chromosomes,
        "output_file": str(combined_output),
        "num_candidates": int(
            sum(int(summary["num_candidates"]) for summary in summaries)
        ),
        "num_independent_snps": int(len(combined)),
        "clump_r2": float(args.clump_r2),
        "clump_window_kb": int(args.clump_kb),
        "clump_populations": args.clump_populations,
        "priority_column": args.priority_column,
        "priority_order": args.priority_order,
        "retention_rule": (
            "retain pairs only when raw r^2 is strictly below clump_r2 "
            "within the physical window"
        ),
        "combined_rank_definition": (
            "CLUMP_RANK is a unique cross-chromosome rank of CLUMP_PRIORITY; "
            "CHROMOSOME_CLUMP_RANK retains the original within-chromosome rank"
        ),
        "chromosome_results": [
            {
                "chromosome": summary["chromosome"],
                "num_candidates": summary["num_candidates"],
                "num_independent_snps": summary["num_independent_snps"],
                "output_file": summary["output_file"],
            }
            for summary in summaries
        ],
        "cross_chromosome_ld_assumed_zero": True,
    }
    combined_summary_path = derived_summary_path(combined_output)
    atomic_write_json(combined_summary, combined_summary_path)
    print(
        f"Combined independent-SNP list: {len(combined):,} rows -> "
        f"{combined_output}"
    )
    print(f"Combined QC: {combined_summary_path}")


def main() -> int:
    args = parse_args()
    ld_dir = args.ld_dir.expanduser().resolve()
    if not ld_dir.is_dir():
        raise FileNotFoundError(f"LD directory does not exist: {ld_dir}")
    output_root = (
        args.out_dir.expanduser().resolve() if args.out_dir is not None else None
    )
    schema_paths = {
        "EUR": (
            args.eur_schema.expanduser().resolve()
            if args.eur_schema is not None
            else (ld_dir / "wb_hapmap/ld_manifest.txt").resolve()
        ),
        "EAS": (
            args.eas_schema.expanduser().resolve()
            if args.eas_schema is not None
            else (ld_dir / "eas_hapmap/ld_manifest.txt").resolve()
        ),
    }
    analysis_groups = (
        [args.clump_populations]
        if args.clump_populations is not None
        else [["EUR"], ["EAS"], ["EUR", "EAS"]]
    )
    if args.combined_output is not None and len(analysis_groups) != 1:
        raise ValueError(
            "--combined-output can only be used with explicit --clump-populations"
        )

    priority_override = args.priority_column
    analyses: list[dict[str, object]] = []
    for populations in analysis_groups:
        args.clump_populations = populations
        args.priority_column = priority_override or (
            "L2_EAS" if populations == ["EAS"] else "L2_EUR"
        )
        for population in populations:
            if not schema_paths[population].is_file():
                raise FileNotFoundError(
                    f"{population} LD manifest does not exist: "
                    f"{schema_paths[population]}"
                )
        jobs = resolve_jobs(
            ld_dir, output_root, args.chromosomes, populations
        )
        combined_output = (
            args.combined_output.expanduser().resolve()
            if args.combined_output is not None
            else default_combined_output(
                ld_dir, output_root, args.chromosomes, populations
            )
        )
        already_complete = preflight(jobs, combined_output, args.overwrite)
        analyses.append(
            {
                "populations": populations,
                "priority_column": args.priority_column,
                "jobs": jobs,
                "combined_output": combined_output,
                "already_complete": already_complete,
            }
        )

    print(f"Chromosomes: {', '.join(args.chromosomes)}")
    print(
        f"Independence rule: raw r^2 < {args.clump_r2:g} within "
        f"+/-{args.clump_kb:,} kb"
    )
    for analysis in analyses:
        populations = analysis["populations"]
        args.clump_populations = populations
        args.priority_column = analysis["priority_column"]
        jobs = analysis["jobs"]
        combined_output = analysis["combined_output"]
        if analysis["already_complete"]:
            print(
                f"Using existing {','.join(populations)} list: "
                f"{combined_output}"
            )
            continue
        print(
            f"\nClumping panels={','.join(populations)}; "
            f"priority={args.priority_column} {args.priority_order}"
        )
        outputs: list[pd.DataFrame] = []
        summaries: list[dict[str, object]] = []
        for job in jobs:
            output, summary = process_chromosome(job, schema_paths, args)
            outputs.append(output)
            summaries.append(summary)
        write_combined_output(outputs, summaries, combined_output, args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
