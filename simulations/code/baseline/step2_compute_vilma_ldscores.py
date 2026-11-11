#!/usr/bin/env python3
"""Compute corrected ViLMA LD scores for one or more chromosomes.

Each chromosome is processed independently from the allele-aligned common-SNP
file produced by ``step1_extract_vilma_snps_bychr.py``. The script computes
within-EUR, within-EAS, and cross-population LD scores, retaining uncorrected
scores for auditing. By default chromosomes 2 and 19 are processed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


DEFAULT_LD_DIR = Path(
    "simulations/data/vilma_ld"
)
DEFAULT_EUR_REFERENCE_N = 10_000
DEFAULT_EAS_REFERENCE_N = 1_154
VAR_COLUMNS = ["ID", "CHROM", "BP", "CM", "A1", "A2"]
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


@dataclass
class SchemaBlock:
    population: str
    index: int
    var_path: Path
    ld_path: Path
    total_variants: int
    selected_positions: np.ndarray
    global_indices: np.ndarray
    signs: np.ndarray
    local_by_global: dict[int, int]


@dataclass
class PreparedSchema:
    population: str
    schema_path: Path
    blocks: dict[int, SchemaBlock]
    membership: np.ndarray
    orientation_counts: dict[str, int]


@dataclass(frozen=True)
class ChromosomeJob:
    chromosome: str
    common_path: Path
    output_path: Path
    summary_path: Path


def file_sha256(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def schema_entries(schema_path: Path):
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


def allele_complement(allele: str) -> str:
    return "".join(COMPLEMENT.get(base, "N") for base in allele.upper())


def allele_orientation(
    canonical_a1: str,
    canonical_a2: str,
    block_a1: str,
    block_a2: str,
) -> tuple[str, int]:
    canonical = (canonical_a1.upper(), canonical_a2.upper())
    block = (block_a1.upper(), block_a2.upper())
    if block == canonical:
        return "same", 1
    if block == canonical[::-1]:
        return "swapped", -1

    complemented = (allele_complement(block[0]), allele_complement(block[1]))
    if complemented == canonical:
        return "complement", 1
    if complemented == canonical[::-1]:
        return "complement_swapped", -1
    raise ValueError(
        f"Alleles cannot be aligned: canonical={canonical}, block={block}"
    )


def load_common_snps(path: Path, expected_snps: int = 0) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Common-SNP file does not exist: {path}")
    common = pd.read_csv(path, sep="\t", low_memory=False)
    required = {"ID", "CHROM", "BP", "A1", "A2"}
    missing = sorted(required - set(common.columns))
    if missing:
        raise ValueError(f"Common-SNP file is missing columns: {missing}")
    if common["ID"].isna().any():
        raise ValueError("Common-SNP file contains missing IDs")
    duplicated = common["ID"].duplicated(keep=False)
    if duplicated.any():
        examples = ", ".join(common.loc[duplicated, "ID"].head(5).astype(str))
        raise ValueError(
            f"Common-SNP file contains {duplicated.sum()} duplicated ID rows; "
            f"examples: {examples}"
        )
    if expected_snps and len(common) != expected_snps:
        raise ValueError(
            f"Expected {expected_snps:,} SNPs but found {len(common):,} in {path}"
        )
    common = common.copy()
    common["ID"] = common["ID"].astype(str)
    common["A1"] = common["A1"].astype(str).str.upper()
    common["A2"] = common["A2"].astype(str).str.upper()
    return common


def prepare_schema(
    schema_path: Path,
    common: pd.DataFrame,
    population: str,
) -> PreparedSchema:
    schema_path = schema_path.expanduser().resolve()
    if not schema_path.is_file():
        raise FileNotFoundError(f"{population} schema does not exist: {schema_path}")

    id_to_global = {snp_id: index for index, snp_id in enumerate(common["ID"])}
    membership = np.full(len(common), -1, dtype=np.int32)
    blocks: dict[int, SchemaBlock] = {}
    orientation_counts: Counter[str] = Counter()

    for block_index, var_path, ld_path in schema_entries(schema_path):
        variants = pd.read_csv(
            var_path,
            sep=r"\s+",
            header=None,
            names=VAR_COLUMNS,
            dtype={"ID": "string", "A1": "string", "A2": "string"},
        )
        selected_positions = np.array(
            [
                position
                for position, snp_id in enumerate(variants["ID"].astype(str))
                if snp_id in id_to_global
            ],
            dtype=np.int64,
        )
        if selected_positions.size == 0:
            continue

        selected = variants.iloc[selected_positions]
        global_indices = np.array(
            [id_to_global[str(snp_id)] for snp_id in selected["ID"]],
            dtype=np.int64,
        )
        duplicate_assignment = membership[global_indices] >= 0
        if duplicate_assignment.any():
            duplicated_ids = common.iloc[global_indices[duplicate_assignment]]["ID"]
            examples = ", ".join(duplicated_ids.head(5).astype(str))
            raise ValueError(
                f"{population} schema assigns SNPs to multiple blocks; "
                f"examples: {examples}"
            )

        signs = np.empty(len(global_indices), dtype=np.float64)
        for local_index, (global_index, block_a1, block_a2) in enumerate(
            zip(global_indices, selected["A1"], selected["A2"])
        ):
            try:
                status, sign = allele_orientation(
                    str(common.iloc[global_index]["A1"]),
                    str(common.iloc[global_index]["A2"]),
                    str(block_a1),
                    str(block_a2),
                )
            except ValueError as exc:
                snp_id = common.iloc[global_index]["ID"]
                raise ValueError(
                    f"{population} allele mismatch for {snp_id}: {exc}"
                ) from exc
            signs[local_index] = sign
            orientation_counts[status] += 1

        membership[global_indices] = block_index
        local_by_global = {
            int(global_index): local_index
            for local_index, global_index in enumerate(global_indices)
        }
        blocks[block_index] = SchemaBlock(
            population=population,
            index=block_index,
            var_path=var_path,
            ld_path=ld_path,
            total_variants=len(variants),
            selected_positions=selected_positions,
            global_indices=global_indices,
            signs=signs,
            local_by_global=local_by_global,
        )

    missing_indices = np.flatnonzero(membership < 0)
    if missing_indices.size:
        examples = ", ".join(
            common.iloc[missing_indices]["ID"].head(5).astype(str)
        )
        raise ValueError(
            f"{population} schema is missing {missing_indices.size:,} common "
            f"SNPs; examples: {examples}"
        )

    expected_block_column = f"{population}_BLOCK"
    if expected_block_column in common.columns:
        expected = pd.to_numeric(
            common[expected_block_column], errors="coerce"
        ).to_numpy(dtype=np.float64)
        mismatch = ~np.isfinite(expected)
        finite = ~mismatch
        mismatch[finite] = expected[finite] != membership[finite]
        if mismatch.any():
            examples = ", ".join(common.loc[mismatch, "ID"].head(5).astype(str))
            raise ValueError(
                f"{population} block assignments disagree with "
                f"{expected_block_column} for {mismatch.sum():,} SNPs; "
                f"examples: {examples}"
            )

    chromosomes = ",".join(
        sorted(common["CHROM"].astype(str).drop_duplicates().tolist())
    )
    print(
        f"Prepared {population}: {len(blocks):,} blocks for chromosome(s) "
        f"{chromosomes}, {len(common):,} SNPs, "
        f"orientations={dict(orientation_counts)}"
    )
    return PreparedSchema(
        population=population,
        schema_path=schema_path,
        blocks=blocks,
        membership=membership,
        orientation_counts=dict(orientation_counts),
    )


def load_aligned_block(block: SchemaBlock) -> np.ndarray:
    stored = np.load(block.ld_path, mmap_mode="r")
    positions = block.selected_positions

    if stored.ndim == 0:
        if block.total_variants != 1:
            raise ValueError(
                f"Scalar LD matrix but {block.total_variants} variants in "
                f"{block.var_path}"
            )
        matrix = np.array([[float(stored)]], dtype=np.float64)
    elif stored.shape[0] == stored.shape[1]:
        if stored.shape[0] != block.total_variants:
            raise ValueError(
                f"LD/variant size mismatch for {block.ld_path}: "
                f"matrix={stored.shape}, variants={block.total_variants}"
            )
        matrix = np.asarray(stored[np.ix_(positions, positions)], dtype=np.float64)
    elif stored.shape[0] > stored.shape[1]:
        num_variants = stored.shape[0] - 1
        if num_variants != block.total_variants:
            raise ValueError(
                f"Low-rank LD/variant size mismatch for {block.ld_path}: "
                f"matrix={stored.shape}, variants={block.total_variants}"
            )
        eigenvectors = np.asarray(stored[:-1, :][positions, :], dtype=np.float64)
        eigenvalues = np.asarray(stored[-1, :], dtype=np.float64)
        matrix = (eigenvectors * eigenvalues) @ eigenvectors.T
    else:
        raise ValueError(f"Unsupported ViLMA LD matrix shape: {stored.shape}")

    matrix *= np.outer(block.signs, block.signs)
    matrix = 0.5 * (matrix + matrix.T)
    if matrix.shape != (len(positions), len(positions)):
        raise AssertionError(f"Unexpected aligned matrix shape for {block.ld_path}")
    if not np.isfinite(matrix).all():
        raise ValueError(f"Non-finite values in LD matrix: {block.ld_path}")
    return matrix


def adjusted_r2(correlation: np.ndarray, reference_n: int) -> np.ndarray:
    """Return the finite-reference approximately unbiased LDSC r-squared."""
    squared = np.square(np.asarray(correlation, dtype=np.float64))
    return squared - (1.0 - squared) / (reference_n - 2.0)


def adjusted_cross_product(
    correlation_1: np.ndarray,
    correlation_2: np.ndarray,
    reference_n_1: int,
    reference_n_2: int,
) -> np.ndarray:
    """Return Popcorn's finite-reference correction for cross-population LD."""
    product = np.asarray(correlation_1, dtype=np.float64) * np.asarray(
        correlation_2, dtype=np.float64
    )
    denominator = np.sqrt(
        (reference_n_1 - 2.0) * (reference_n_2 - 2.0)
    )
    correction = (1.0 - np.abs(product)) / denominator
    return product - np.sign(product) * correction


def compute_within_population_scores(
    prepared: PreparedSchema,
    num_snps: int,
    progress_every: int,
    reference_n: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    corrected_scores = np.full(num_snps, np.nan, dtype=np.float64)
    raw_scores = np.full(num_snps, np.nan, dtype=np.float64)
    diagonals = np.full(num_snps, np.nan, dtype=np.float64)
    blocks = sorted(prepared.blocks.values(), key=lambda block: block.index)

    for count, block in enumerate(blocks, start=1):
        matrix = load_aligned_block(block)
        raw_scores[block.global_indices] = np.einsum(
            "ij,ij->i", matrix, matrix, optimize=True
        )
        corrected_scores[block.global_indices] = np.sum(
            adjusted_r2(matrix, reference_n), axis=1
        )
        diagonals[block.global_indices] = np.diag(matrix)
        if progress_every and (count % progress_every == 0 or count == len(blocks)):
            print(
                f"{prepared.population} LD scores: processed "
                f"{count:,}/{len(blocks):,} blocks"
            )

    if (
        not np.isfinite(corrected_scores).all()
        or not np.isfinite(raw_scores).all()
        or not np.isfinite(diagonals).all()
    ):
        raise ValueError(f"Incomplete {prepared.population} LD scores")
    return corrected_scores, raw_scores, diagonals


def compute_cross_population_scores(
    eur: PreparedSchema,
    eas: PreparedSchema,
    num_snps: int,
    progress_every: int,
    eur_reference_n: int,
    eas_reference_n: int,
) -> tuple[np.ndarray, np.ndarray, int]:
    groups: dict[tuple[int, int], list[int]] = defaultdict(list)
    for global_index, pair in enumerate(zip(eur.membership, eas.membership)):
        groups[(int(pair[0]), int(pair[1]))].append(global_index)

    all_blocks = {"EUR": eur.blocks, "EAS": eas.blocks}

    @lru_cache(maxsize=6)
    def cached_matrix(population: str, block_index: int) -> np.ndarray:
        return load_aligned_block(all_blocks[population][block_index])

    corrected_cross_scores = np.full(num_snps, np.nan, dtype=np.float64)
    raw_cross_scores = np.full(num_snps, np.nan, dtype=np.float64)
    sorted_groups = sorted(groups.items())
    for count, ((eur_block_index, eas_block_index), indices) in enumerate(
        sorted_groups, start=1
    ):
        global_indices = np.asarray(indices, dtype=np.int64)
        eur_block = eur.blocks[eur_block_index]
        eas_block = eas.blocks[eas_block_index]
        eur_local = np.array(
            [eur_block.local_by_global[int(index)] for index in global_indices],
            dtype=np.int64,
        )
        eas_local = np.array(
            [eas_block.local_by_global[int(index)] for index in global_indices],
            dtype=np.int64,
        )

        eur_matrix = cached_matrix("EUR", eur_block_index)
        eas_matrix = cached_matrix("EAS", eas_block_index)
        eur_overlap = eur_matrix[np.ix_(eur_local, eur_local)]
        eas_overlap = eas_matrix[np.ix_(eas_local, eas_local)]
        raw_cross_scores[global_indices] = np.einsum(
            "ij,ij->i", eur_overlap, eas_overlap, optimize=True
        )
        corrected_cross_scores[global_indices] = np.sum(
            adjusted_cross_product(
                eur_overlap,
                eas_overlap,
                eur_reference_n,
                eas_reference_n,
            ),
            axis=1,
        )

        if progress_every and (
            count % progress_every == 0 or count == len(sorted_groups)
        ):
            print(
                "Cross-population LD scores: processed "
                f"{count:,}/{len(sorted_groups):,} overlapping block pairs"
            )

    cached_matrix.cache_clear()
    if (
        not np.isfinite(corrected_cross_scores).all()
        or not np.isfinite(raw_cross_scores).all()
    ):
        raise ValueError("Incomplete cross-population LD scores")
    return corrected_cross_scores, raw_cross_scores, len(sorted_groups)


def numeric_summary(values: np.ndarray) -> dict[str, float]:
    return {
        "min": float(np.min(values)),
        "q01": float(np.quantile(values, 0.01)),
        "median": float(np.median(values)),
        "mean": float(np.mean(values)),
        "q99": float(np.quantile(values, 0.99)),
        "max": float(np.max(values)),
    }


def derived_summary_path(output: Path) -> Path:
    text = str(output)
    for suffix in (".tsv.gz", ".tsv"):
        if text.endswith(suffix):
            return Path(text[: -len(suffix)] + ".summary.json")
    return Path(text + ".summary.json")


def normalized_chromosome(value: object) -> str:
    chromosome = re.sub(r"^CHR", "", str(value).strip().upper())
    if re.fullmatch(r"0*[0-9]+(?:\.0+)?", chromosome):
        return str(int(float(chromosome)))
    return chromosome


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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compute corrected EUR, EAS, and cross-population ViLMA LD scores "
            "for one or more chromosomes."
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
            "Optional output root. Chromosome-specific chrN_snps directories "
            "are created below it; default: write beside each input file"
        ),
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
        "--eur-reference-n",
        type=int,
        default=DEFAULT_EUR_REFERENCE_N,
        help="EUR/White British LD-reference sample size (default: 10000)",
    )
    parser.add_argument(
        "--eas-reference-n",
        type=int,
        default=DEFAULT_EAS_REFERENCE_N,
        help="EAS LD-reference sample size (default: 1154)",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=20,
        help="Print progress after this many blocks; use 0 to disable",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace existing LD-score and summary files",
    )
    args = parser.parse_args()
    try:
        args.chromosomes = normalized_chromosomes(args.chromosomes)
    except ValueError as exc:
        parser.error(str(exc))
    if args.eur_reference_n <= 2:
        parser.error("--eur-reference-n must be greater than 2")
    if args.eas_reference_n <= 2:
        parser.error("--eas-reference-n must be greater than 2")
    if args.progress_every < 0:
        parser.error("--progress-every cannot be negative")
    return args


def resolve_jobs(
    ld_dir: Path,
    output_root: Path | None,
    chromosomes: list[str],
) -> list[ChromosomeJob]:
    jobs: list[ChromosomeJob] = []
    for chromosome in chromosomes:
        input_dir = ld_dir / f"chr{chromosome}_snps"
        common_path = (
            input_dir / f"eas_eur_chr{chromosome}_common_snps.tsv"
        ).resolve()
        output_dir = (
            input_dir
            if output_root is None
            else output_root / f"chr{chromosome}_snps"
        )
        output_path = (
            output_dir / f"eas_eur_chr{chromosome}_ldscores.tsv.gz"
        ).resolve()
        jobs.append(
            ChromosomeJob(
                chromosome=chromosome,
                common_path=common_path,
                output_path=output_path,
                summary_path=derived_summary_path(output_path),
            )
        )
    return jobs


def preflight_jobs(jobs: list[ChromosomeJob], overwrite: bool) -> None:
    missing = [job.common_path for job in jobs if not job.common_path.is_file()]
    if missing:
        missing_text = "\n  ".join(str(path) for path in missing)
        raise FileNotFoundError(f"Missing Step 1 common-SNP file(s):\n  {missing_text}")

    if overwrite:
        return
    existing = [
        path
        for job in jobs
        for path in (job.output_path, job.summary_path)
        if path.exists()
    ]
    if existing:
        existing_text = "\n  ".join(str(path) for path in existing)
        raise FileExistsError(
            f"Refusing to overwrite existing output(s):\n  {existing_text}\n"
            "Use --overwrite to replace them."
        )


def validate_common_chromosome(
    chromosome: str, common_path: Path, common
) -> None:
    observed = {
        normalized_chromosome(value)
        for value in common["CHROM"].dropna().unique()
    }
    if observed != {chromosome}:
        raise ValueError(
            f"{common_path} contains chromosome labels {sorted(observed)}; "
            f"expected only chromosome {chromosome}"
        )


def build_output(common, scores: dict[str, np.ndarray]):
    output_columns = [
        column
        for column in ("ID", "CHROM", "BP", "A1", "A2", "AMBIGUOUS")
        if column in common.columns
    ]
    output = common.loc[:, output_columns].copy()
    for column, values in scores.items():
        output[column] = values
    output["L2.pop1"] = scores["L2_EUR"]
    output["L2.pop2"] = scores["L2_EAS"]
    output["L12"] = scores["L12_EUR_EAS"]
    return output


def write_result(output, summary: dict[str, object], job: ChromosomeJob) -> None:
    job.output_path.parent.mkdir(parents=True, exist_ok=True)
    job.summary_path.parent.mkdir(parents=True, exist_ok=True)
    output_tmp = job.output_path.with_name(job.output_path.name + ".tmp")
    summary_tmp = job.summary_path.with_name(job.summary_path.name + ".tmp")
    output.to_csv(
        output_tmp,
        sep="\t",
        index=False,
        float_format="%.12g",
        compression="gzip",
    )
    summary_tmp.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    output_tmp.replace(job.output_path)
    summary_tmp.replace(job.summary_path)


def compute_chromosome(
    job: ChromosomeJob,
    requested_chromosomes: list[str],
    eur_schema_path: Path,
    eas_schema_path: Path,
    eur_schema_sha256: str,
    eas_schema_sha256: str,
    args: argparse.Namespace,
) -> dict[str, object]:
    common = load_common_snps(job.common_path, expected_snps=0)
    if common.empty:
        raise ValueError(f"Step 1 common-SNP file is empty: {job.common_path}")
    validate_common_chromosome(job.chromosome, job.common_path, common)
    print(
        f"Chromosome {job.chromosome}: loaded {len(common):,} common SNPs "
        f"from {job.common_path}"
    )

    eur = prepare_schema(eur_schema_path, common, "EUR")
    eas = prepare_schema(eas_schema_path, common, "EAS")
    l2_eur, l2_eur_raw, diag_eur = compute_within_population_scores(
        eur,
        len(common),
        args.progress_every,
        args.eur_reference_n,
    )
    l2_eas, l2_eas_raw, diag_eas = compute_within_population_scores(
        eas,
        len(common),
        args.progress_every,
        args.eas_reference_n,
    )
    l12, l12_raw, num_block_pairs = compute_cross_population_scores(
        eur,
        eas,
        len(common),
        args.progress_every,
        args.eur_reference_n,
        args.eas_reference_n,
    )

    if (l2_eur_raw < -1e-12).any() or (l2_eas_raw < -1e-12).any():
        raise ValueError(
            f"Chromosome {job.chromosome}: negative raw within-population "
            "LD score detected"
        )
    raw_cauchy_bound = np.sqrt(
        np.maximum(l2_eur_raw, 0.0) * np.maximum(l2_eas_raw, 0.0)
    )
    max_raw_cauchy_excess = float(
        np.max(np.abs(l12_raw) - raw_cauchy_bound)
    )
    if max_raw_cauchy_excess > 1e-8:
        raise ValueError(
            f"Chromosome {job.chromosome}: raw cross-population LD score "
            "violates the Cauchy-Schwarz bound; maximum excess="
            f"{max_raw_cauchy_excess:.6g}"
        )

    scores = {
        "L2_EUR": l2_eur,
        "L2_EAS": l2_eas,
        "L12_EUR_EAS": l12,
        "L2_EUR_RAW": l2_eur_raw,
        "L2_EAS_RAW": l2_eas_raw,
        "L12_EUR_EAS_RAW": l12_raw,
    }
    output = build_output(common, scores)
    summary: dict[str, object] = {
        "chromosome": job.chromosome,
        "requested_chromosomes": requested_chromosomes,
        "num_snps": int(len(common)),
        "common_snp_file": str(job.common_path),
        "common_snp_file_sha256": file_sha256(job.common_path),
        "output_file": str(job.output_path),
        "population_mapping": {"pop1": "EUR/White British", "pop2": "EAS"},
        "ld_reference_sample_sizes": {
            "EUR_White_British": int(args.eur_reference_n),
            "EAS": int(args.eas_reference_n),
        },
        "eur_schema": str(eur_schema_path),
        "eas_schema": str(eas_schema_path),
        "eur_schema_sha256": eur_schema_sha256,
        "eas_schema_sha256": eas_schema_sha256,
        "num_eur_blocks": int(len(eur.blocks)),
        "num_eas_blocks": int(len(eas.blocks)),
        "num_overlapping_block_pairs": int(num_block_pairs),
        "eur_orientation_counts": eur.orientation_counts,
        "eas_orientation_counts": eas.orientation_counts,
        "score_definitions": {
            "L2.pop1": (
                "sum over chromosome-specific common SNPs of "
                "R_EUR[j,k]^2 - (1 - R_EUR[j,k]^2) / (N_EUR - 2)"
            ),
            "L2.pop2": (
                "sum over chromosome-specific common SNPs of "
                "R_EAS[j,k]^2 - (1 - R_EAS[j,k]^2) / (N_EAS - 2)"
            ),
            "L12": (
                "sum over chromosome-specific common SNPs of "
                "c - sign(c) * (1 - abs(c)) / "
                "sqrt((N_EUR - 2) * (N_EAS - 2)), where "
                "c = R_EUR[j,k] * R_EAS[j,k]"
            ),
            "L2_EUR_RAW": "sum of R_EUR[j,k]^2",
            "L2_EAS_RAW": "sum of R_EAS[j,k]^2",
            "L12_EUR_EAS_RAW": "sum of R_EUR[j,k] * R_EAS[j,k]",
            "includes_self_correlation": True,
            "allele_reference": f"A1/A2 in {job.common_path.name}",
        },
        "finite_reference_bias_correction": {
            "within_population": (
                "LDSC approximately unbiased r-squared estimator: "
                "r2 - (1 - r2) / (N_ref - 2)"
            ),
            "cross_population": (
                "Popcorn correction for independent reference panels: "
                "c - sign(c) * (1 - abs(c)) / "
                "sqrt((N_ref1 - 2) * (N_ref2 - 2))"
            ),
            "negative_corrected_scores_are_allowed": True,
        },
        "ld_matrix_handling": (
            "ViLMA square matrices and stacked low-rank decompositions are "
            "supported; scores use reconstructed matrices without diagonal "
            "renormalization. Finite-reference formulas treat reconstructed "
            "entries as estimated correlations."
        ),
        "L2_EUR_summary": numeric_summary(l2_eur),
        "L2_EAS_summary": numeric_summary(l2_eas),
        "L12_summary": numeric_summary(l12),
        "L2_EUR_raw_summary": numeric_summary(l2_eur_raw),
        "L2_EAS_raw_summary": numeric_summary(l2_eas_raw),
        "L12_raw_summary": numeric_summary(l12_raw),
        "fraction_L2_EUR_corrected_below_zero": float(np.mean(l2_eur < 0)),
        "fraction_L2_EAS_corrected_below_zero": float(np.mean(l2_eas < 0)),
        "fraction_L2_EUR_corrected_below_one": float(np.mean(l2_eur < 1)),
        "fraction_L2_EAS_corrected_below_one": float(np.mean(l2_eas < 1)),
        "EUR_diagonal_summary": numeric_summary(diag_eur),
        "EAS_diagonal_summary": numeric_summary(diag_eas),
        "max_raw_cauchy_schwarz_excess": max_raw_cauchy_excess,
    }
    write_result(output, summary, job)
    print(f"Wrote chromosome {job.chromosome} LD scores: {job.output_path}")
    print(f"Wrote chromosome {job.chromosome} QC: {job.summary_path}")
    print(
        f"Chromosome {job.chromosome} corrected means: "
        f"EUR={np.mean(l2_eur):.6g}, EAS={np.mean(l2_eas):.6g}, "
        f"cross={np.mean(l12):.6g}"
    )
    return {
        "chromosome": job.chromosome,
        "num_snps": int(len(common)),
        "output_file": str(job.output_path),
        "summary_file": str(job.summary_path),
    }


def main() -> int:
    args = parse_args()
    ld_dir = args.ld_dir.expanduser().resolve()
    if not ld_dir.is_dir():
        raise FileNotFoundError(f"LD directory does not exist: {ld_dir}")
    output_root = (
        args.out_dir.expanduser().resolve() if args.out_dir is not None else None
    )
    eur_schema_path = (
        args.eur_schema.expanduser().resolve()
        if args.eur_schema is not None
        else (ld_dir / "wb_hapmap/ld_manifest.txt").resolve()
    )
    eas_schema_path = (
        args.eas_schema.expanduser().resolve()
        if args.eas_schema is not None
        else (ld_dir / "eas_hapmap/ld_manifest.txt").resolve()
    )
    for population, schema_path in (
        ("EUR/White British", eur_schema_path),
        ("EAS", eas_schema_path),
    ):
        if not schema_path.is_file():
            raise FileNotFoundError(
                f"{population} LD manifest does not exist: {schema_path}"
            )

    jobs = resolve_jobs(ld_dir, output_root, args.chromosomes)
    preflight_jobs(jobs, args.overwrite)
    eur_schema_sha256 = file_sha256(eur_schema_path)
    eas_schema_sha256 = file_sha256(eas_schema_path)
    print(f"Chromosomes: {', '.join(args.chromosomes)}")
    print(f"EUR/White British schema: {eur_schema_path}")
    print(f"EAS schema: {eas_schema_path}")
    print(
        "LD reference sample sizes: "
        f"EUR={args.eur_reference_n:,}, EAS={args.eas_reference_n:,}"
    )

    results = [
        compute_chromosome(
            job,
            args.chromosomes,
            eur_schema_path,
            eas_schema_path,
            eur_schema_sha256,
            eas_schema_sha256,
            args,
        )
        for job in jobs
    ]
    total_snps = sum(int(result["num_snps"]) for result in results)
    print(
        f"Completed {len(results)} chromosome(s), {total_snps:,} total "
        "chromosome-specific common SNP rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
