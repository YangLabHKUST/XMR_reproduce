#!/usr/bin/env python3
"""Generate XMR true effects after randomly selecting one fixed SNP architecture.

Core candidates meet both EUR and EAS LD-score 45th-75th percentile intervals.
Select 50 shared, 5 target-only and 5 auxiliary-only cores by default in random order:
shared+auxiliary cores are independent in EUR and shared+target cores in EAS
(10 Mb, raw r^2 < 0.001). Expand once in either population at 1 Mb/r^2 >= 0.4.
Every core and expanded neighbor has foreground AND background effects.

Select background-only SNPs comprising 70% of each chromosome, independent of ALL
expanded foreground in both populations (1 Mb, raw r^2 < 0.001).
Foreground is added to this 70% when counting M_bg. Background-only SNPs may
have mutual LD. SNP identities/classes stay fixed
across replicates and settings; only Gaussian effects and GWAS errors are redrawn.
Background h2 is 0.2 for exposures and 0.3 for outcomes; active gamma
variance is 10 times the exposure-background per-SNP variance.
"""

from __future__ import annotations

import os
import argparse
import hashlib
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from step2_compute_vilma_ldscores import (
    PreparedSchema,
    load_aligned_block,
    prepare_schema,
)


XMR_SOURCE = (
    "https://github.com/YangLabHKUST/XMR_reproduce/blob/main/"
    "sim/code/simulation_reproduce_plot.ipynb"
)
DEFAULT_LD_DIR = Path(
    "simulations/data/vilma_ld"
)
DEFAULT_OUT_DIR = Path(
    "simulations/data/"
    "reduced_confounding/true_effects"
)
DEFAULT_EUR_Y_H2 = 0.3
DEFAULT_EUR_X_H2 = 0.2
DEFAULT_EAS_X_H2 = 0.2
DEFAULT_EAS_Y_H2 = 0.3
DEFAULT_BACKGROUND_ONLY_FRACTION = 0.70
DEFAULT_N_TARGET_ONLY = 5
DEFAULT_N_AUX_ONLY = 5
DEFAULT_N_SHARED = 50
DEFAULT_N_INDEPENDENT_FOREGROUND = DEFAULT_N_SHARED
DEFAULT_GAMMA_VARIANCE_RATIO = 10.0
DEFAULT_ALPHA_VARIANCE_RATIO = 1.0
DEFAULT_N_AUX_VALUES = (40_000, 60_000)
DEFAULT_N_TARGET = 10_000
ARCHITECTURE_VERSION = "fixed_mid_ldscore_bgonly_v4_master_seed"
CORE_POPULATIONS = {
    "shared": ("EUR", "EAS"), "target_only": ("EAS",), "aux_only": ("EUR",),
}

# Order: V1, U1, U2, V2. The U1/U2 entry is replaced by --rhos.
# EUR and EAS background X/Y correlations are fixed at 0.4.
R_OMEGA_TEMPLATE = np.array(
    [
        [1.0, 0.4, 0.1, 0.7],
        [0.4, 1.0, 0.5, 0.3],
        [0.1, 0.5, 1.0, 0.4],
        [0.7, 0.3, 0.4, 1.0],
    ],
    dtype=np.float64,
)

# Sampling-error covariance in summary-statistic order Y1, X1, X2, Y2.
# All EUR (Y1, X1) versus EAS (X2, Y2) intercepts are zero.
C_MATRIX = np.array(
    [
        [1.10, 0.15, 0.0, 0.0],
        [0.15, 1.25, 0.0, 0.0],
        [0.0, 0.0, 1.10, 0.15],
        [0.0, 0.0, 0.15, 1.05],
    ],
    dtype=np.float64,
)
# EAS order: X2, Y2; configurable independently of foreground rho.
C_MATRIX[2, 2] = float(os.environ.get("XMR_EAS_C11", "1.10"))
C_MATRIX[2, 3] = C_MATRIX[3, 2] = float(os.environ.get("XMR_EAS_C12", "0.15"))
C_MATRIX[3, 3] = float(os.environ.get("XMR_EAS_C22", "1.05"))
if not np.isfinite(C_MATRIX).all() or np.linalg.eigvalsh(C_MATRIX).min() <= 0:
    raise ValueError("Configured C must be finite and positive definite")



@dataclass(frozen=True)
class Architecture:
    foreground_indices: np.ndarray
    independent_foreground_indices: np.ndarray
    additional_foreground_indices: np.ndarray
    target_only_indices: np.ndarray
    aux_only_indices: np.ndarray
    shared_indices: np.ndarray
    background_indices: np.ndarray
    additional_background_indices: np.ndarray
    xmr_order: np.ndarray
    has_background: np.ndarray
    is_foreground: np.ndarray
    is_core_foreground: np.ndarray
    ld_neighbor_of_core: np.ndarray
    iv_aux: np.ndarray
    iv_target: np.ndarray
    effect_class: np.ndarray
    foreground_source: np.ndarray
    background_by_chromosome: dict[str, dict[str, int | float]]


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


def parse_float_list(raw: str | Iterable[object], option: str) -> list[float]:
    values: list[float] = []
    fields = [raw] if isinstance(raw, str) else raw
    for field in fields:
        for item in str(field).split(","):
            item = item.strip()
            if not item:
                continue
            try:
                value = float(item)
            except ValueError as exc:
                raise argparse.ArgumentTypeError(
                    f"{option} contains a non-numeric value: {item!r}"
                ) from exc
            if not math.isfinite(value):
                raise argparse.ArgumentTypeError(
                    f"{option} values must be finite: {item!r}"
                )
            if value not in values:
                values.append(value)
    if not values:
        raise argparse.ArgumentTypeError(f"{option} must contain a value")
    return values


def parse_positive_int_list(
    raw: str | Iterable[object], option: str
) -> list[int]:
    numeric = parse_float_list(raw, option)
    values: list[int] = []
    for value in numeric:
        integer = int(value)
        if value != integer or integer <= 0:
            raise argparse.ArgumentTypeError(
                f"{option} values must be positive integers: {value:g}"
            )
        if integer not in values:
            values.append(integer)
    return values


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate XMR true effects with one randomly selected, fixed chr2/chr19 architecture."
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
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument(
        "--eur-ldscore-column",
        default="L2_EUR",
        help="EUR LD-score column retained in outputs",
    )
    parser.add_argument(
        "--n-independent-foreground",
        type=int,
        default=None,
        help="Shared core count (legacy alias); defaults to --n-shared and must equal it",
    )
    parser.add_argument(
        "--n-target-only", type=int, default=DEFAULT_N_TARGET_ONLY,
        help="Number of target-only core SNPs, before LD expansion",
    )
    parser.add_argument("--n-aux-only", type=int, default=DEFAULT_N_AUX_ONLY,
                        help="Number of auxiliary-only core SNPs, before LD expansion")
    parser.add_argument("--n-shared", type=int, default=DEFAULT_N_SHARED,
                        help="Number of shared core SNPs, before LD expansion")
    parser.add_argument(
        "--background-only-fraction",
        type=float,
        default=DEFAULT_BACKGROUND_ONLY_FRACTION,
        help="Background-only fraction of all SNPs per chromosome; foreground is additional (default: 0.70)",
    )
    parser.add_argument(
        "--h2",
        dest="h2_common",
        type=float,
        default=None,
        help="Legacy common h2 override for all four components",
    )
    parser.add_argument(
        "--h2-outcome",
        type=float,
        default=None,
        help="Legacy common outcome h2 override for V1 and V2",
    )
    parser.add_argument(
        "--h2-eur-y",
        type=float,
        default=None,
        help="EUR Y background h2; also scales ALPHA1 (default: 0.3)",
    )
    parser.add_argument(
        "--h2-eur-x",
        type=float,
        default=None,
        help="EUR X background h2; also scales GAMMA1 (default: 0.2)",
    )
    parser.add_argument(
        "--h2-eas-x",
        type=float,
        default=None,
        help="EAS X background h2; also scales GAMMA2 (default: 0.2)",
    )
    parser.add_argument(
        "--h2-eas-y",
        type=float,
        default=None,
        help="EAS Y background h2; also scales ALPHA2 (default: 0.3)",
    )
    parser.add_argument(
        "--foreground-gamma-variance-ratio",
        type=float,
        default=DEFAULT_GAMMA_VARIANCE_RATIO,
        help="Foreground gamma/background variance ratio (default: 10)",
    )
    parser.add_argument(
        "--foreground-alpha-variance-ratio",
        type=float,
        default=DEFAULT_ALPHA_VARIANCE_RATIO,
    )
    parser.add_argument("--ldscore-quantiles", type=float, nargs=2, default=[0.45, 0.75])
    parser.add_argument("--core-kb", type=int, default=10000)
    parser.add_argument("--core-r2", type=float, default=0.001)
    parser.add_argument("--ld-kb", type=int, default=1000, help="Foreground expansion window")
    parser.add_argument("--ld-r2", type=float, default=0.4, help="Foreground expansion raw r^2")
    parser.add_argument("--background-exclusion-kb", type=int, default=1000)
    parser.add_argument("--background-exclusion-r2", type=float, default=0.001)
    parser.add_argument(
        "--n-aux",
        nargs="+",
        default=[str(value) for value in DEFAULT_N_AUX_VALUES],
        help=(
            "One or more auxiliary sample sizes, separated by spaces or "
            "commas. Each value creates a setting directory for Step 5; true "
            "effects do not depend on sample size (default: 40000 60000)."
        ),
    )
    parser.add_argument(
        "--n-target",
        nargs="+",
        default=[str(DEFAULT_N_TARGET)],
        help=(
            "One or more target sample sizes, separated by spaces or commas. "
            "Each value creates a setting directory for Step 5; true effects "
            "do not depend on sample size (default: 10000)."
        ),
    )
    parser.add_argument(
        "--target",
        default="EAS",
        help="Target population label used in output names (default: EAS)",
    )
    parser.add_argument(
        "--betas",
        nargs="+",
        default=["0", "0.2"],
        help=(
            "One or more beta values separated by spaces or commas; "
            "beta1=beta2 (default: 0 0.2)"
        ),
    )
    parser.add_argument(
        "--rhos",
        nargs="+",
        default=["0.4", "0.8"],
        help=(
            "One or more correlations separated by spaces or commas for both "
            "U1/U2 and GAMMA1/GAMMA2 (default: 0.4 0.8)"
        ),
    )
    parser.add_argument("--replicates", type=int, default=50)
    parser.add_argument("--replicate-start", type=int, default=1)
    parser.add_argument(
        "--seed", type=int, default=414,
        help="Master seed; component- and replicate-specific streams are derived from it",
    )
    parser.add_argument(
        "--paired-effect-streams", action="store_true",
        help=("Pair background and foreground Gaussian draws by original SNP index "
              "across background fractions; use a fresh output root"),
    )
    parser.add_argument(
        "--expected-total-snps",
        type=int,
        default=0,
        help="Fail unless the combined Step 2 inputs contain this many SNPs",
    )
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.n_independent_foreground is None:
        args.n_independent_foreground = args.n_shared

    try:
        args.chromosomes = normalized_chromosomes(args.chromosomes)
        args.betas = parse_float_list(args.betas, "--betas")
        args.rhos = parse_float_list(args.rhos, "--rhos")
        args.n_aux = parse_positive_int_list(args.n_aux, "--n-aux")
        args.n_target = parse_positive_int_list(args.n_target, "--n-target")
    except (ValueError, argparse.ArgumentTypeError) as exc:
        parser.error(str(exc))
    common_h2 = args.h2_common
    outcome_h2 = (
        args.h2_outcome if args.h2_outcome is not None else common_h2
    )
    args.h2_eur_y = (
        args.h2_eur_y
        if args.h2_eur_y is not None
        else (outcome_h2 if outcome_h2 is not None else DEFAULT_EUR_Y_H2)
    )
    args.h2_eur_x = (
        args.h2_eur_x
        if args.h2_eur_x is not None
        else (common_h2 if common_h2 is not None else DEFAULT_EUR_X_H2)
    )
    args.h2_eas_x = (
        args.h2_eas_x
        if args.h2_eas_x is not None
        else (common_h2 if common_h2 is not None else DEFAULT_EAS_X_H2)
    )
    args.h2_eas_y = (
        args.h2_eas_y
        if args.h2_eas_y is not None
        else (outcome_h2 if outcome_h2 is not None else DEFAULT_EAS_Y_H2)
    )
    count_options = {
        "--n-independent-foreground": args.n_independent_foreground,
        "--n-target-only": args.n_target_only,
        "--n-aux-only": args.n_aux_only,
        "--n-shared": args.n_shared,
    }
    for option, value in count_options.items():
        if value < 0:
            parser.error(f"{option} cannot be negative")
    n_foreground = args.n_target_only + args.n_aux_only + args.n_shared
    if n_foreground < 1:
        parser.error("The foreground class counts must sum to at least 1")
    if args.n_independent_foreground > n_foreground:
        parser.error(
            "--n-independent-foreground cannot exceed the total foreground count"
        )
    if args.n_independent_foreground > args.n_shared:
        parser.error(
            "--n-independent-foreground cannot exceed --n-shared because "
            "every anchor is fixed as shared"
        )
    if args.n_independent_foreground != args.n_shared:
        parser.error(
            "All shared cores are anchors; set "
            "--n-independent-foreground equal to --n-shared"
        )
    if not 0 < args.background_only_fraction <= 1:
        parser.error("--background-only-fraction must be in (0, 1]")
    h2_options = {
        "--h2-eur-y": args.h2_eur_y,
        "--h2-eur-x": args.h2_eur_x,
        "--h2-eas-x": args.h2_eas_x,
        "--h2-eas-y": args.h2_eas_y,
    }
    for option, value in h2_options.items():
        if not 0 < value < 1:
            parser.error(f"{option} must be in (0, 1)")
    if args.foreground_gamma_variance_ratio <= 0:
        parser.error("--foreground-gamma-variance-ratio must be positive")
    if args.foreground_alpha_variance_ratio <= 0:
        parser.error("--foreground-alpha-variance-ratio must be positive")
    for name in ("core", "ld", "background_exclusion"):
        if getattr(args, name + "_kb") < 1 or not 0 < getattr(args, name + "_r2") <= 1:
            parser.error(f"{name}: require kb >= 1 and 0 < r2 <= 1")
    q_low, q_high = args.ldscore_quantiles
    if not 0 <= q_low < q_high <= 1:
        parser.error("Require 0 <= LD-score quantile LOW < HIGH <= 1")
    if args.replicates < 1:
        parser.error("--replicates must be at least 1")
    if args.replicate_start < 1:
        parser.error("--replicate-start must be at least 1")
    if not 0 <= args.seed <= 2_147_483_646:
        parser.error("--seed must be between 0 and 2147483646")
    if args.expected_total_snps < 0:
        parser.error("--expected-total-snps cannot be negative")
    if any(abs(rho) >= 1 for rho in args.rhos):
        parser.error("Every --rhos value must satisfy -1 < rho < 1")
    for rho in args.rhos:
        minimum = float(np.linalg.eigvalsh(background_correlation(rho)).min())
        if minimum <= 0:
            parser.error(
                f"rho={rho:g} makes R_Omega non-positive-definite "
                f"(minimum eigenvalue {minimum:.6g})"
            )
    return args


def file_sha256(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def ids_sha256(snps: pd.DataFrame, indices: np.ndarray) -> str:
    text = "\n".join(snps.iloc[indices]["ID"].astype(str))
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def safe_population_label(value: str) -> str:
    label = re.sub(r"[^A-Za-z0-9]+", "_", value.strip().upper()).strip("_")
    if not label:
        raise ValueError("Target population label is empty")
    return label


def float_label(value: float) -> str:
    return (
        format(value, ".12g")
        .replace("-", "m")
        .replace(".", "p")
        .replace("+", "")
    )


def chromosome_label(chromosomes: list[str]) -> str:
    return "_".join(f"chr{chromosome}" for chromosome in chromosomes)


def load_step2_snps(
    ld_dir: Path,
    chromosomes: list[str],
    ldscore_column: str,
    expected_total_snps: int,
) -> tuple[pd.DataFrame, list[dict[str, object]]]:
    frames: list[pd.DataFrame] = []
    inputs: list[dict[str, object]] = []
    required = {
        "ID",
        "CHROM",
        "BP",
        "A1",
        "A2",
        ldscore_column,
        "L2_EUR",
        "L2_EAS",
    }
    for chromosome_index, chromosome in enumerate(chromosomes):
        path = (
            ld_dir
            / f"chr{chromosome}_snps"
            / f"eas_eur_chr{chromosome}_ldscores.tsv.gz"
        ).resolve()
        if not path.is_file():
            raise FileNotFoundError(f"Missing Step 2 LD-score file: {path}")
        frame = pd.read_csv(path, sep="\t", low_memory=False)
        missing = sorted(required - set(frame.columns))
        if missing:
            raise ValueError(f"{path} is missing required columns: {missing}")
        observed = {
            normalized_chromosome(value)
            for value in frame["CHROM"].dropna().unique()
        }
        if observed != {chromosome}:
            raise ValueError(
                f"{path} contains chromosomes {sorted(observed)}; "
                f"expected only {chromosome}"
            )
        frame = frame.copy()
        frame["ID"] = frame["ID"].astype(str)
        frame["CHROM"] = chromosome
        frame["A1"] = frame["A1"].astype(str).str.upper()
        frame["A2"] = frame["A2"].astype(str).str.upper()
        frame["BP"] = pd.to_numeric(frame["BP"], errors="coerce")
        frame[ldscore_column] = pd.to_numeric(
            frame[ldscore_column], errors="coerce"
        )
        if not np.isfinite(frame["BP"].to_numpy(dtype=np.float64)).all():
            raise ValueError(f"Non-finite BP values in {path}")
        if not np.isfinite(
            frame[ldscore_column].to_numpy(dtype=np.float64)
        ).all():
            raise ValueError(f"Non-finite {ldscore_column} values in {path}")
        frame["BP"] = frame["BP"].astype(np.int64)
        frame["CHROMOSOME_INPUT_ORDER"] = chromosome_index
        frames.append(frame)
        inputs.append(
            {
                "chromosome": chromosome,
                "path": str(path),
                "sha256": file_sha256(path),
                "num_snps": int(len(frame)),
            }
        )
    snps = pd.concat(frames, ignore_index=True, sort=False)
    if snps["ID"].isna().any():
        raise ValueError("Step 2 inputs contain missing SNP IDs")
    duplicated = snps["ID"].duplicated(keep=False)
    if duplicated.any():
        examples = ", ".join(snps.loc[duplicated, "ID"].head(5))
        raise ValueError(
            f"Step 2 inputs contain {duplicated.sum():,} duplicate ID rows; "
            f"examples: {examples}"
        )
    if expected_total_snps and len(snps) != expected_total_snps:
        raise ValueError(
            f"Expected {expected_total_snps:,} total SNPs but found "
            f"{len(snps):,}"
        )
    return snps, inputs


def core_candidates(
    snps: pd.DataFrame, quantiles: list[float]
) -> tuple[np.ndarray, dict[str, dict[str, float | str]]]:
    """Intersect inclusive EUR/EAS quantile intervals, calculated over all input SNPs."""
    candidate_mask = np.ones(len(snps), dtype=bool)
    bounds = {}
    for population in ("EUR", "EAS"):
        column = f"L2_{population}"
        scores = pd.to_numeric(snps[column], errors="coerce").to_numpy(dtype=float)
        if not np.isfinite(scores).all():
            raise ValueError(f"Non-finite values in {column}")
        lower, upper = np.quantile(scores, quantiles)
        candidate_mask &= (scores >= lower) & (scores <= upper)
        bounds[population] = {"column": column, "lower": float(lower), "upper": float(upper)}
    return np.flatnonzero(candidate_mask), bounds


def select_random_cores(
    snps: pd.DataFrame,
    candidate_indices: np.ndarray,
    sizes: dict[str, int],
    schemas: dict[str, PreparedSchema],
    rng: np.random.Generator,
    window_bp: int,
    r2_threshold: float,
) -> dict[str, np.ndarray]:
    """Random-order greedy selection; never relax thresholds or retry failed draws."""
    cores = {}
    used = np.zeros(len(snps), dtype=bool)
    blocked = {population: np.zeros(len(snps), dtype=bool) for population in schemas}
    for label, populations in CORE_POPULATIONS.items():
        chosen = []
        if sizes[label]:
            for index in rng.permutation(candidate_indices[~used[candidate_indices]]):
                if any(blocked[population][index] for population in populations):
                    continue
                chosen.append(int(index))
                used[index] = True
                for population in populations:
                    neighborhood, _ = anchor_ld_neighborhood(
                        np.array([index]), snps, {population: schemas[population]},
                        window_bp, r2_threshold,
                    )
                    blocked[population] |= neighborhood
                if len(chosen) == sizes[label]:
                    break
        if len(chosen) < sizes[label]:
            raise ValueError(
                f"{label}: need {sizes[label]}, retained {len(chosen)} after random-order LD selection"
            )
        cores[label] = np.asarray(chosen, dtype=np.int64)
    return cores


def anchor_ld_neighborhood(
    anchor_indices: np.ndarray,
    snps: pd.DataFrame,
    schemas: dict[str, PreparedSchema],
    window_bp: int,
    r2_threshold: float,
) -> tuple[np.ndarray, dict[str, object]]:
    """Return the union of EUR/EAS SNPs linked to at least one anchor."""
    anchor_indices = np.asarray(anchor_indices, dtype=np.int64)
    if np.unique(anchor_indices).size != len(anchor_indices):
        raise ValueError("Anchor indices must be unique")
    bp = snps["BP"].to_numpy(dtype=np.int64)
    chromosome = snps["CHROM"].astype(str).to_numpy()
    union = np.zeros(len(snps), dtype=bool)
    by_population: dict[str, dict[str, int | float]] = {}

    for population, prepared in schemas.items():
        population_mask = np.zeros(len(snps), dtype=bool)
        memberships = prepared.membership[anchor_indices]
        maximum_linked_r2 = 0.0
        linked_pairs = 0
        for block_index in np.unique(memberships):
            block = prepared.blocks[int(block_index)]
            block_anchors = anchor_indices[memberships == block_index]
            matrix = load_aligned_block(block)
            block_globals = block.global_indices
            block_positions = bp[block_globals]
            nonanchor = ~np.isin(block_globals, anchor_indices)
            linked_in_block = np.zeros(len(block_globals), dtype=bool)
            for anchor in block_anchors:
                anchor_local = block.local_by_global[int(anchor)]
                in_window = (
                    (chromosome[block_globals] == chromosome[anchor])
                    & (np.abs(block_positions - bp[anchor]) <= window_bp)
                )
                r2 = np.square(matrix[anchor_local, :])
                linked = in_window & (r2 >= r2_threshold)
                linked_in_block |= linked
                linked_nonanchor = linked & nonanchor
                linked_pairs += int(np.count_nonzero(linked_nonanchor))
                if linked_nonanchor.any():
                    maximum_linked_r2 = max(
                        maximum_linked_r2, float(np.max(r2[linked_nonanchor]))
                    )
            population_mask[block_globals[linked_in_block]] = True
        union |= population_mask
        nonanchor_mask = population_mask.copy()
        nonanchor_mask[anchor_indices] = False
        by_population[population] = {
            "num_snps_including_anchors": int(np.count_nonzero(population_mask)),
            "num_nonanchor_snps_in_neighborhood": int(np.count_nonzero(nonanchor_mask)),
            "num_anchor_to_nonanchor_linked_pairs": linked_pairs,
            "maximum_linked_pair_r2": maximum_linked_r2,
        }

    union[anchor_indices] = True
    nonanchor_union = union.copy()
    nonanchor_union[anchor_indices] = False
    summary: dict[str, object] = {
        "populations": list(schemas),
        "r2_threshold": float(r2_threshold),
        "window_kb": int(window_bp // 1_000),
        "window_bp": int(window_bp),
        "neighborhood_rule": (
            "include a nonanchor SNP if raw r^2 meets or exceeds the threshold "
            "with any supplied anchor in any supplied population within the physical window"
        ),
        "by_population": by_population,
        "num_anchors": int(len(anchor_indices)),
        "num_nonanchor_snps_in_union": int(
            np.count_nonzero(nonanchor_union)
        ),
        "num_snps_including_anchors_in_union": int(
            np.count_nonzero(union)
        ),
        "cross_block_correlations_assumed_zero": True,
    }
    return union, summary


def rounded_fraction_count(total: int, fraction: float) -> int:
    return int(math.floor(total * fraction + 0.5))


def expand_foreground_from_cores(
    snps: pd.DataFrame,
    core_groups: dict[str, np.ndarray],
    schemas: dict[str, PreparedSchema],
    window_bp: int,
    r2_threshold: float,
) -> tuple[dict[str, np.ndarray], np.ndarray, dict[str, object]]:
    """One-hop expansion; core labels are fixed and neighbor labels are merged."""
    cores = np.concatenate(list(core_groups.values()))
    if np.unique(cores).size != len(cores):
        raise ValueError("Core foreground groups overlap")
    masks = {}
    details = {}
    for label, indices in core_groups.items():
        masks[label], details[label] = anchor_ld_neighborhood(
            indices, snps, schemas, window_bp, r2_threshold
        )
        masks[label][cores] = False
    active_aux = masks["shared"] | masks["aux_only"]
    active_target = masks["shared"] | masks["target_only"]
    neighbor_classes = {
        "shared": active_aux & active_target,
        "target_only": active_target & ~active_aux,
        "aux_only": active_aux & ~active_target,
    }
    groups = {
        label: np.concatenate([core_groups[label], np.flatnonzero(mask)])
        for label, mask in neighbor_classes.items()
    }
    region = active_aux | active_target
    region[cores] = True
    metadata = {
        "rule": "one_hop_from_core_no_recursion",
        "neighbor_class_rule": "union_of_linked_core_active_populations_core_classes_fixed",
        "effect_draws_across_snps": "independent_same_per_snp_foreground_variance",
        "populations": list(schemas),
        "r2_threshold": float(r2_threshold),
        "window_kb": window_bp // 1000,
        "num_core": int(len(cores)),
        "num_ld_neighbors": int(region.sum() - len(cores)),
        "num_foreground": int(region.sum()),
        "expanded_class_counts": {label: len(indices) for label, indices in groups.items()},
        "by_core_class": details,
        "cross_block_correlations_assumed_zero": True,
    }
    return groups, region, metadata


def select_background(
    snps: pd.DataFrame,
    foreground_indices: np.ndarray,
    foreground_ld_mask: np.ndarray,
    background_exclusion_mask: np.ndarray,
    chromosomes: list[str],
    background_only_fraction: float,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray, dict[str, dict[str, int | float]]]:
    """Draw the background-only quota outside all expanded foreground LD neighborhoods."""
    foreground_indices = np.asarray(foreground_indices, dtype=np.int64)
    if np.unique(foreground_indices).size != len(foreground_indices):
        raise ValueError("Foreground indices must be unique")

    chromosome_values = snps["CHROM"].astype(str).to_numpy()
    foreground_mask = np.zeros(len(snps), dtype=bool)
    foreground_mask[foreground_indices] = True
    foreground_ld_mask = np.asarray(foreground_ld_mask, dtype=bool)
    if foreground_ld_mask.shape != (len(snps),):
        raise ValueError("Foreground LD-neighborhood mask must match the SNP list")
    if not foreground_ld_mask[foreground_indices].all():
        raise ValueError("LD-neighborhood mask must include all foreground SNPs")
    background_exclusion_mask = np.asarray(background_exclusion_mask, dtype=bool)
    if background_exclusion_mask.shape != (len(snps),) or not background_exclusion_mask[foreground_ld_mask].all():
        raise ValueError("Background exclusion must cover the expanded foreground set")
    selected_by_chromosome: list[np.ndarray] = []
    summary: dict[str, dict[str, int | float]] = {}

    for chromosome in chromosomes:
        chromosome_indices = np.flatnonzero(chromosome_values == chromosome)
        chromosome_foreground = chromosome_indices[
            foreground_mask[chromosome_indices]
        ]
        target_count = rounded_fraction_count(
            len(chromosome_indices), background_only_fraction
        )
        chromosome_neighbors = chromosome_indices[
            foreground_ld_mask[chromosome_indices] & ~foreground_mask[chromosome_indices]
        ]
        num_ld_neighbors = len(chromosome_neighbors)
        candidates = chromosome_indices[~background_exclusion_mask[chromosome_indices]]
        needed = target_count
        if needed > len(candidates):
            raise ValueError(
                f"Chromosome {chromosome}: need {needed:,} background-only SNPs "
                f"for the {background_only_fraction:.1%} background-only quota, but only {len(candidates):,} "
                "are LD-independent of ALL expanded foreground SNPs in EUR/EAS. "
                f"Shortfall={needed - len(candidates):,}. No LD thresholds were relaxed."
            )
        random_background = rng.permutation(candidates)[:needed]
        selected = np.concatenate([
            chromosome_foreground, chromosome_neighbors, random_background
        ])
        selected_by_chromosome.append(selected)
        summary[chromosome] = {
            "num_snps": int(len(chromosome_indices)),
            "num_core_forced_into_background": int(
                len(chromosome_foreground)
            ),
            "num_random_nonforeground_background": int(len(random_background)),
            "num_eligible_isolated_background": int(len(candidates)),
            "num_isolated_background": int(len(random_background)),
            "num_excluded_buffer_snps": int(np.count_nonzero(
                background_exclusion_mask[chromosome_indices] & ~foreground_ld_mask[chromosome_indices]
            )),
            "num_ld_neighbors": int(num_ld_neighbors),
            "num_background": int(len(selected)),
            "num_union": int(len(selected)),
            "realized_union_fraction": float(len(selected) / len(chromosome_indices)),
            "requested_background_only_count": int(target_count),
            "realized_background_fraction": float(len(selected) / len(chromosome_indices)),
            "realized_background_only_fraction": float(len(random_background) / len(chromosome_indices)),
        }

    background_indices = np.concatenate(selected_by_chromosome)
    if np.unique(background_indices).size != len(background_indices):
        raise AssertionError("Background selection produced duplicate SNPs")
    background_mask = np.zeros(len(snps), dtype=bool)
    background_mask[background_indices] = True
    if not background_mask[foreground_indices].all():
        raise AssertionError("Every foreground SNP must belong to background")
    if not background_mask[foreground_ld_mask].all():
        raise AssertionError("Every foreground LD neighbor must have a background effect")
    if np.any(background_mask & ~foreground_ld_mask & background_exclusion_mask):
        raise AssertionError("Background-only overlaps the expanded foreground LD neighborhood")
    return background_indices, background_mask, summary


def build_architecture(
    snps: pd.DataFrame,
    core_groups: dict[str, np.ndarray],
    foreground_groups: dict[str, np.ndarray],
    background_indices: np.ndarray,
    background_by_chromosome: dict[str, dict[str, int | float]],
    rng: np.random.Generator,
) -> Architecture:
    num_snps = len(snps)
    shared_indices = foreground_groups["shared"]
    target_only = foreground_groups["target_only"]
    aux_only = foreground_groups["aux_only"]
    core_indices = np.concatenate(list(core_groups.values()))
    core_mask = np.zeros(num_snps, dtype=bool)
    core_mask[core_indices] = True
    foreground_order = np.concatenate([shared_indices, target_only, aux_only])
    additional_indices = foreground_order[~np.isin(foreground_order, core_groups["shared"])]
    if np.unique(foreground_order).size != len(foreground_order):
        raise ValueError("Foreground classes overlap by SNP index")

    foreground_mask = np.zeros(num_snps, dtype=bool)
    foreground_mask[foreground_order] = True
    background_indices = np.asarray(background_indices, dtype=np.int64)
    if np.unique(background_indices).size != len(background_indices):
        raise ValueError("Background indices must be unique")
    background_mask = np.zeros(num_snps, dtype=bool)
    background_mask[background_indices] = True
    if not background_mask[foreground_order].all():
        raise AssertionError("Every foreground SNP must have a background effect")
    if not foreground_mask[core_indices].all():
        raise AssertionError("Every core SNP must belong to expanded foreground")
    ld_neighbors = foreground_mask & ~core_mask

    nonforeground_background = background_indices[
        ~foreground_mask[background_indices]
    ]
    background_order = np.concatenate(
        [foreground_order, nonforeground_background]
    )
    if np.unique(background_order).size != len(background_indices):
        raise AssertionError("Background architecture contains duplicate SNPs")

    chromosome_values = snps["CHROM"].astype(str).to_numpy()
    background_summary = {
        chromosome: dict(values)
        for chromosome, values in background_by_chromosome.items()
    }
    for chromosome, values in background_summary.items():
        chromosome_mask = chromosome_values == chromosome
        values["num_core"] = int(np.count_nonzero(chromosome_mask & core_mask))
        values["num_foreground"] = int(
            np.count_nonzero(chromosome_mask & foreground_mask)
        )
        values["num_nonforeground_background"] = int(
            np.count_nonzero(
                chromosome_mask & background_mask & ~foreground_mask
            )
        )

    null_indices = np.flatnonzero(~background_mask)
    null_order = rng.permutation(null_indices)
    full_order = np.concatenate(
        [foreground_order, nonforeground_background, null_order]
    )
    if np.unique(full_order).size != num_snps:
        raise AssertionError("XMR order must contain every SNP exactly once")
    xmr_order = np.empty(num_snps, dtype=np.int64)
    xmr_order[full_order] = np.arange(1, num_snps + 1, dtype=np.int64)

    has_background = background_mask.astype(np.int8)
    is_foreground = foreground_mask.astype(np.int8)
    iv_aux = np.zeros(num_snps, dtype=np.int8)
    iv_target = np.zeros(num_snps, dtype=np.int8)
    iv_target[target_only] = 1
    iv_aux[aux_only] = 1
    iv_aux[shared_indices] = 1
    iv_target[shared_indices] = 1
    effect_class = np.full(num_snps, "null_effect", dtype="<U16")
    effect_class[background_indices] = "background"
    effect_class[target_only] = "target_only"
    effect_class[aux_only] = "aux_only"
    effect_class[shared_indices] = "shared"
    foreground_source = np.full(num_snps, "not_foreground", dtype="<U32")
    foreground_source[ld_neighbors] = "core_ld_neighbor"
    foreground_source[core_groups["shared"]] = "random_shared_core"
    foreground_source[core_groups["target_only"]] = "random_target_only_core"
    foreground_source[core_groups["aux_only"]] = "random_aux_only_core"
    return Architecture(
        foreground_indices=foreground_order,
        independent_foreground_indices=core_groups["shared"],
        additional_foreground_indices=additional_indices,
        target_only_indices=target_only,
        aux_only_indices=aux_only,
        shared_indices=shared_indices,
        background_indices=background_order,
        additional_background_indices=nonforeground_background,
        xmr_order=xmr_order,
        has_background=has_background,
        is_foreground=is_foreground,
        is_core_foreground=core_mask.astype(np.int8),
        ld_neighbor_of_core=ld_neighbors.astype(np.int8),
        iv_aux=iv_aux,
        iv_target=iv_target,
        effect_class=effect_class,
        foreground_source=foreground_source,
        background_by_chromosome=background_summary,
    )


def background_correlation(rho: float) -> np.ndarray:
    correlation = R_OMEGA_TEMPLATE.copy()
    background_rho = rho
    eas_corr = float(os.environ.get("XMR_EAS_BACKGROUND_CORR", "0.4"))
    correlation[1, 2] = correlation[2, 1] = background_rho
    correlation[2, 3] = correlation[3, 2] = eas_corr
    correlation[1, 3] = correlation[3, 1] = float(os.environ.get("XMR_BG_EURX_EASY_CORR", "0.3"))
    correlation[0, 2] = correlation[2, 0] = float(os.environ.get("XMR_BG_EURY_EASX_CORR", "0.1"))
    if not np.isfinite(correlation).all() or np.linalg.eigvalsh(correlation).min() <= 0:
        raise ValueError("Configured background correlation must be finite and positive definite")
    return correlation


def background_covariance(
    rho: float, omega_variances: np.ndarray
) -> np.ndarray:
    scales = np.sqrt(np.asarray(omega_variances, dtype=np.float64))
    if scales.shape != (4,) or not np.isfinite(scales).all():
        raise ValueError("omega_variances must contain four finite values")
    return np.outer(scales, scales) * background_correlation(rho)


def foreground_covariance(
    rho: float,
    omega_variances: np.ndarray,
    gamma_variance_ratio: float,
    alpha_variance_ratio: float,
) -> np.ndarray:
    omega_variances = np.asarray(omega_variances, dtype=np.float64)
    variances = omega_variances * np.array(
        [
            alpha_variance_ratio,
            gamma_variance_ratio,
            gamma_variance_ratio,
            alpha_variance_ratio,
        ],
        dtype=np.float64,
    )
    correlation = np.eye(4, dtype=np.float64)
    correlation[1, 2] = rho
    correlation[2, 1] = rho
    scales = np.sqrt(variances)
    return np.outer(scales, scales) * correlation


def covariance_sqrt(matrix: np.ndarray, name: str) -> np.ndarray:
    eigenvalues = np.linalg.eigvalsh(matrix)
    if eigenvalues.min() <= 0:
        raise ValueError(
            f"{name} is not positive definite; minimum eigenvalue "
            f"{eigenvalues.min():.6g}"
        )
    return np.linalg.cholesky(matrix)


def replicate_latents(
    architecture: Architecture, rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    return (
        rng.standard_normal((len(architecture.background_indices), 4)),
        rng.standard_normal((len(architecture.foreground_indices), 4)),
    )


def effect_seed(master_seed: int, replicate: int) -> np.random.SeedSequence:
    """Derive an order-independent true-effect stream for one replicate."""
    return np.random.SeedSequence([master_seed, 2, replicate])


def paired_replicate_latents(
    architecture: Architecture, master_seed: int, replicate: int, num_snps: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Common SNPs keep their draws when the background-only fraction changes.

    Draw in original Step 2 SNP order before subsetting. Separate background and
    foreground streams prevent a changed M_bg from shifting foreground draws.
    The default production RNG remains in replicate_latents for compatibility.
    """
    background_rng = np.random.default_rng(
        np.random.SeedSequence([master_seed, 2, replicate, 0])
    )
    foreground_rng = np.random.default_rng(
        np.random.SeedSequence([master_seed, 2, replicate, 1])
    )
    return (
        background_rng.standard_normal((num_snps, 4))[architecture.background_indices],
        foreground_rng.standard_normal((num_snps, 4))[architecture.foreground_indices],
    )


def generate_background(
    architecture: Architecture,
    latent: np.ndarray,
    rho: float,
    omega_variances: np.ndarray,
    num_snps: int,
) -> np.ndarray:
    values = latent @ covariance_sqrt(
        background_covariance(rho, omega_variances), "Omega"
    ).T
    background = np.zeros((num_snps, 4), dtype=np.float64)
    background[architecture.background_indices] = values
    return background


def generate_foreground(
    architecture: Architecture,
    latent: np.ndarray,
    rho: float,
    omega_variances: np.ndarray,
    gamma_variance_ratio: float,
    alpha_variance_ratio: float,
    num_snps: int,
) -> np.ndarray:
    sigma_f = foreground_covariance(
        rho,
        omega_variances,
        gamma_variance_ratio,
        alpha_variance_ratio,
    )
    values = latent @ covariance_sqrt(sigma_f, "Sigma_f").T
    foreground = np.zeros((num_snps, 4), dtype=np.float64)
    foreground[architecture.foreground_indices] = values
    foreground[architecture.target_only_indices, 1] = 0.0
    foreground[architecture.aux_only_indices, 2] = 0.0
    return foreground


def metadata_columns(snps: pd.DataFrame) -> list[str]:
    preferred = ["ID", "CHROM", "BP", "A1", "A2", "AMBIGUOUS"]
    score_columns = [
        column
        for column in snps.columns
        if column.startswith("L2") or column.startswith("L12")
    ]
    columns = [column for column in preferred if column in snps.columns]
    columns.extend(column for column in score_columns if column not in columns)
    return columns


def build_output(
    snps: pd.DataFrame,
    architecture: Architecture,
    background: np.ndarray,
    foreground: np.ndarray,
    beta2: float,
) -> pd.DataFrame:
    output = snps.loc[:, metadata_columns(snps)].copy()
    output["XMR_ORDER"] = architecture.xmr_order
    output["EFFECT_CLASS"] = architecture.effect_class
    output["FOREGROUND_SOURCE"] = architecture.foreground_source
    output["HAS_BACKGROUND"] = architecture.has_background
    output["Z_FOREGROUND"] = architecture.is_foreground
    output["IS_CORE_FOREGROUND"] = architecture.is_core_foreground
    output["LD_NEIGHBOR_OF_CORE"] = architecture.ld_neighbor_of_core
    output["IN_ARCHITECTURE_UNION"] = architecture.has_background
    output["IV_AUX"] = architecture.iv_aux
    output["IV_TARGET"] = architecture.iv_target
    output["V1"] = background[:, 0]
    output["U1"] = background[:, 1]
    output["U2"] = background[:, 2]
    output["V2"] = background[:, 3]
    output["ALPHA1"] = foreground[:, 0]
    output["GAMMA1"] = foreground[:, 1]
    output["GAMMA2"] = foreground[:, 2]
    output["ALPHA2"] = foreground[:, 3]
    output["ALPHA"] = output["ALPHA2"]
    output["TRUE_X1"] = output["U1"] + output["GAMMA1"]
    output["TRUE_X2"] = output["U2"] + output["GAMMA2"]
    beta1 = beta2
    output["TRUE_Y1"] = (
        output["V1"] + output["ALPHA1"] + beta1 * output["GAMMA1"]
    )
    output["TRUE_Y2"] = (
        output["V2"] + output["ALPHA2"] + beta2 * output["GAMMA2"]
    )
    return output


def safe_correlation(x: np.ndarray, y: np.ndarray) -> float | None:
    if len(x) < 2 or np.std(x) == 0 or np.std(y) == 0:
        return None
    return float(np.corrcoef(x, y)[0, 1])


def empirical_qc(
    snps: pd.DataFrame,
    architecture: Architecture,
    background: np.ndarray,
    foreground: np.ndarray,
) -> dict[str, object]:
    background_values = background[architecture.background_indices]
    shared_values = foreground[architecture.shared_indices]
    counts_by_chromosome: dict[str, dict[str, int]] = {}
    chromosome_values = snps["CHROM"].astype(str).to_numpy()
    for chromosome in pd.unique(chromosome_values):
        mask = chromosome_values == chromosome
        counts_by_chromosome[str(chromosome)] = {
            "num_snps": int(np.count_nonzero(mask)),
            "num_background": int(
                np.count_nonzero(mask & architecture.has_background.astype(bool))
            ),
            "num_foreground": int(
                np.count_nonzero(mask & architecture.is_foreground.astype(bool))
            ),
        }
    return {
        "effect_class_counts": {
            str(key): int(value)
            for key, value in pd.Series(architecture.effect_class)
            .value_counts()
            .items()
        },
        "counts_by_chromosome": counts_by_chromosome,
        "all_foreground_nested_in_background": bool(
            np.all(
                architecture.has_background[architecture.foreground_indices] == 1
            )
        ),
        "background_empirical_correlation": np.corrcoef(
            background_values, rowvar=False
        ).tolist(),
        "background_exposure_empirical_correlation": safe_correlation(
            background_values[:, 1], background_values[:, 2]
        ),
        "shared_foreground_gamma_empirical_correlation": safe_correlation(
            shared_values[:, 1], shared_values[:, 2]
        ),
        "sum_squared_true_components": {
            "V1": float(np.square(background[:, 0]).sum()),
            "U1": float(np.square(background[:, 1]).sum()),
            "U2": float(np.square(background[:, 2]).sum()),
            "V2": float(np.square(background[:, 3]).sum()),
            "ALPHA1": float(np.square(foreground[:, 0]).sum()),
            "GAMMA1": float(np.square(foreground[:, 1]).sum()),
            "GAMMA2": float(np.square(foreground[:, 2]).sum()),
            "ALPHA2": float(np.square(foreground[:, 3]).sum()),
        },
    }


def scenario_label(args: argparse.Namespace) -> str:
    return (
        f"h2eury{float_label(args.h2_eur_y)}_eurx{float_label(args.h2_eur_x)}_"
        f"easx{float_label(args.h2_eas_x)}_easy{float_label(args.h2_eas_y)}_"
        f"bgOnly{float_label(100 * args.background_only_fraction)}pct_"
        f"core{args.n_shared}s{args.n_target_only}t{args.n_aux_only}a_"
        f"q{float_label(100 * args.ldscore_quantiles[0])}to{float_label(100 * args.ldscore_quantiles[1])}_"
        f"c{float_label(args.core_kb / 1000)}MbR2lt{float_label(args.core_r2)}_"
        f"e{float_label(args.ld_kb / 1000)}MbR2ge{float_label(args.ld_r2)}_"
        f"b{float_label(args.background_exclusion_kb / 1000)}MbR2lt"
        f"{float_label(args.background_exclusion_r2)}_fixedIV"
    )


def setting_id(beta: float, rho: float, n_aux: int, n_target: int) -> str:
    return (
        f"beta{float_label(beta)}_rho{float_label(rho)}_"
        f"naux{n_aux}_ntarget{n_target}"
    )


def output_path_for(
    out_dir: Path,
    target: str,
    chromosomes: list[str],
    rho: float,
    beta: float,
    n_aux: int,
    n_target: int,
    gamma_ratio: float,
    scenario: str,
    replicate: int,
) -> Path:
    stem = (
        f"xmr_truth_WB_{target}_{chromosome_label(chromosomes)}_"
        f"rho{float_label(rho)}_beta{float_label(beta)}_"
        f"gvar{float_label(gamma_ratio)}_{scenario}_rep{replicate:03d}"
    )
    return out_dir / setting_id(beta, rho, n_aux, n_target) / f"{stem}.tsv.gz"


def metadata_path_for(output_path: Path) -> Path:
    text = str(output_path)
    if not text.endswith(".tsv.gz"):
        raise ValueError(f"Output path must end in .tsv.gz: {output_path}")
    return Path(text[: -len(".tsv.gz")] + ".json")


def atomic_write_table(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    frame.to_csv(
        temporary,
        sep="\t",
        index=False,
        compression="gzip",
        float_format="%.17g",
    )
    temporary.replace(path)


def atomic_write_json(data: dict[str, object], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary.replace(path)


def output_pair_is_complete(
    data_path: Path, metadata_path: Path, overwrite: bool
) -> bool:
    if overwrite:
        return False
    if data_path.exists() and metadata_path.exists():
        return True
    if data_path.exists() or metadata_path.exists():
        raise FileExistsError(
            "Incomplete checkpoint; remove it or use --overwrite:\n  "
            f"{data_path}\n  {metadata_path}"
        )
    return False


def architecture_table(
    snps: pd.DataFrame,
    architecture: Architecture,
) -> pd.DataFrame:
    frame = snps.loc[:, metadata_columns(snps)].copy()
    frame["XMR_ORDER"] = architecture.xmr_order
    frame["EFFECT_CLASS"] = architecture.effect_class
    frame["FOREGROUND_SOURCE"] = architecture.foreground_source
    frame["HAS_BACKGROUND"] = architecture.has_background
    frame["Z_FOREGROUND"] = architecture.is_foreground
    frame["IS_CORE_FOREGROUND"] = architecture.is_core_foreground
    frame["LD_NEIGHBOR_OF_CORE"] = architecture.ld_neighbor_of_core
    frame["IN_ARCHITECTURE_UNION"] = architecture.has_background
    frame["IV_AUX"] = architecture.iv_aux
    frame["IV_TARGET"] = architecture.iv_target
    return frame


def save_checkpoint(
    frame: pd.DataFrame, metadata: dict[str, object], path: Path, overwrite: bool,
) -> None:
    metadata_path = metadata_path_for(path)
    if output_pair_is_complete(path, metadata_path, overwrite):
        saved = json.loads(metadata_path.read_text(encoding="utf-8"))
        for key in ("design_id", "replicate", "architecture_id", "setting_id"):
            if saved.get(key) != metadata.get(key):
                raise ValueError(f"{path}: checkpoint {key} differs; use a new output directory")
        if saved.get("output_file_sha256") != file_sha256(path):
            raise ValueError(f"{path}: checkpoint data changed; use --overwrite to regenerate")
        print(f"Checked existing: {path}", flush=True)
        return
    atomic_write_table(frame, path)
    metadata = dict(metadata, output_file_sha256=file_sha256(path))
    atomic_write_json(metadata, metadata_path)
    print(f"Wrote: {path}", flush=True)


def main() -> int:
    args = parse_args()
    ld_dir = args.ld_dir.expanduser().resolve()
    if not ld_dir.is_dir():
        raise FileNotFoundError(f"LD directory does not exist: {ld_dir}")
    out_dir = args.out_dir.expanduser().resolve()
    if any("_fixedIV_rep" not in path.name for path in out_dir.rglob("xmr_truth_*.tsv.gz")):
        raise ValueError(
            "This output directory contains an older SNP architecture. Use a new --out-dir "
            "(and a new study --work-root) so old summary/background/method results are not reused."
        )
    target = safe_population_label(args.target)

    # 1. Read existing LD scores, form the mid-LD-score candidate set, and load LD.
    snps, step2_inputs = load_step2_snps(
        ld_dir, args.chromosomes, args.eur_ldscore_column, args.expected_total_snps,
    )
    sizes = {"shared": args.n_shared, "target_only": args.n_target_only, "aux_only": args.n_aux_only}
    candidates, score_bounds = core_candidates(snps, args.ldscore_quantiles)
    if len(candidates) < sum(sizes.values()):
        raise ValueError(f"Only {len(candidates)} SNPs meet both LD-score intervals; need {sum(sizes.values())} cores")
    schema_paths = {
        "EUR": (args.eur_schema or ld_dir / "wb_hapmap/ld_manifest.txt").expanduser().resolve(),
        "EAS": (args.eas_schema or ld_dir / "eas_hapmap/ld_manifest.txt").expanduser().resolve(),
    }
    schemas = {
        population: prepare_schema(path, snps, population)
        for population, path in schema_paths.items()
    }
    ld_metadata = {
        population: {"path": str(path), "sha256": file_sha256(path),
                     "num_blocks_used": len(schemas[population].blocks)}
        for population, path in schema_paths.items()
    }
    h2_background = np.array([args.h2_eur_y, args.h2_eur_x, args.h2_eas_x, args.h2_eas_y])
    component_names = ("V1", "U1", "U2", "V2")
    scenario = scenario_label(args)
    configuration = {
        "architecture_version": ARCHITECTURE_VERSION,
        "scenario": scenario, "seed": args.seed, "target": target,
        "step2_inputs": step2_inputs, "ld_schemas": ld_metadata,
        "gamma_variance_ratio": args.foreground_gamma_variance_ratio,
        "alpha_variance_ratio": args.foreground_alpha_variance_ratio,
        "C": C_MATRIX.tolist(), "R_Omega_template": R_OMEGA_TEMPLATE.tolist(),
    }
    # Do not change existing design hashes when the optional mode is disabled.
    if args.paired_effect_streams:
        configuration["effect_stream_mode"] = "paired_by_original_snp"
    design_id = hashlib.sha256(json.dumps(configuration, sort_keys=True).encode()).hexdigest()
    design_path = out_dir / "simulation_design.json"
    if design_path.is_file():
        previous = json.loads(design_path.read_text(encoding="utf-8"))
        if previous.get("design_id") != design_id:
            raise ValueError(
                f"{design_path}: simulation design changed; use a new --out-dir/work root "
                "to avoid reusing incompatible downstream results"
            )
    atomic_write_json(dict(configuration, design_id=design_id), design_path)
    print(f"Loaded {len(snps):,} SNPs; eligible cores={len(candidates):,}; counts={sizes}", flush=True)
    print(f"Core LD-score bounds: {score_bounds}", flush=True)
    print("SNP architecture is selected once and fixed across all replicates/settings.", flush=True)

    # Core selection matches the first draw of the diagnostic with the same seed.
    core_rng = np.random.default_rng(np.random.SeedSequence([args.seed, 0, 1]))
    background_rng = np.random.default_rng(np.random.SeedSequence([args.seed, 1]))
    # 2. Select cores, expand foreground once, then draw background-only SNPs.
    core_groups = select_random_cores(
        snps, candidates, sizes, schemas, core_rng, args.core_kb * 1000, args.core_r2,
    )
    core_indices = np.concatenate(list(core_groups.values()))
    groups, foreground_mask, expansion = expand_foreground_from_cores(
        snps, core_groups, schemas, args.ld_kb * 1000, args.ld_r2,
    )
    exclusion_mask, exclusion = anchor_ld_neighborhood(
        np.flatnonzero(foreground_mask), snps, schemas,
        args.background_exclusion_kb * 1000, args.background_exclusion_r2,
    )
    background_indices, _, background_by_chr = select_background(
        snps, core_indices, foreground_mask, exclusion_mask,
        args.chromosomes, args.background_only_fraction, background_rng,
    )
    architecture = build_architecture(
        snps, core_groups, groups, background_indices, background_by_chr, background_rng,
    )
    num_foreground = len(architecture.foreground_indices)
    num_background = len(architecture.background_indices)
    omega_variances = h2_background / num_background
    architecture_frame = architecture_table(snps, architecture)
    architecture_id = hashlib.sha256(
        pd.util.hash_pandas_object(architecture_frame, index=False).to_numpy().tobytes()
    ).hexdigest()
    counts = [{
        "n_core": len(core_indices), "n_foreground": num_foreground,
        "n_ld_neighbors": int(architecture.ld_neighbor_of_core.sum()),
        **{f"n_{label}": len(indices) for label, indices in groups.items()},
        "n_background": num_background, "n_background_only": num_background - num_foreground,
        "architecture_id": architecture_id,
    }]

    architecture_path = out_dir / "architectures" / (
        f"{scenario}_seed{args.seed}_architecture.tsv.gz"
    )
    architecture_metadata = {
        "architecture_version": ARCHITECTURE_VERSION,
        "design_id": design_id, "architecture_id": architecture_id,
        "master_seed": args.seed,
        "base_seed": args.seed,
        "rng_strategy": "one_master_seed_with_component_and_replicate_substreams",
        "rng_streams": {
            "core": "default_rng(SeedSequence([base_seed, 0, 1]))",
            "background_and_null_order": "default_rng(SeedSequence([base_seed, 1]))",
            "effects": "default_rng(SeedSequence([base_seed, 2, replicate]))",
        },
        "paired_across_settings_within_replicate": True,
        "chromosomes": args.chromosomes, "num_snps": len(snps),
        "step2_inputs": step2_inputs, "ld_schemas": ld_metadata,
        "core_required_populations": {key: list(value) for key, value in CORE_POPULATIONS.items()},
        "foreground_selection": {
            "num_core": len(core_indices), "num_foreground": num_foreground,
            "core_counts": sizes, "ld_expansion": expansion,
            "selection_rule": "random_order_greedy_LD_selection_without_replacement_no_retries",
            "ldscore_quantiles": args.ldscore_quantiles,
            "ldscore_bounds": score_bounds,
            "ldscore_quantile_reference": "all_input_SNPs_pooled_across_chromosomes_per_population",
            "ldscore_filter": "intersection_of_both_inclusive_intervals_core_only",
            "num_core_candidates": len(candidates),
            "groups_disjoint_by_snp_id": True,
            "independence_threshold_r2": args.core_r2,
            "independence_window_kb": args.core_kb,
            "selected_all_id_sha256": ids_sha256(snps, np.sort(architecture.foreground_indices)),
            "selected_core_id_sha256": {
                key: ids_sha256(snps, value) for key, value in core_groups.items()
            },
            "fixed_across_replicates": True,
        },
        "foreground_class_assignment": {
            **{key: len(value) for key, value in groups.items()},
            **{f"core_{key}": value for key, value in sizes.items()},
            "groups_disjoint_by_snp_id": True,
            "assignment_id_order_sha256": ids_sha256(snps, architecture.foreground_indices),
            "fixed_across_replicates": True,
        },
        "background_selection": {
            "method": "include_all_expanded_foreground_then_random_isolated_background_only",
            "fraction_denominator": "all_common_SNPs_per_chromosome",
            "fraction_numerator": "background_only",
            "M_bg_definition": "expanded_foreground_plus_background_only",
            "ld_neighbors_have_background_effects": True,
            "ld_neighbors_have_foreground_effects": True,
            "isolated_background_foreground_ld_independent": True,
            "background_only_mutual_independence_required": False,
            "foreground_ld_neighborhood": expansion,
            "background_only_ld_exclusion": exclusion,
            "background_isolation_reference": "all_expanded_foreground",
            "requested_background_only_fraction_per_chromosome": args.background_only_fraction,
            "rounding": "background_only_count=floor(M_chr * fraction + 0.5)",
            "by_chromosome": architecture.background_by_chromosome,
            "num_background": num_background,
            "num_ld_neighbors": int(architecture.ld_neighbor_of_core.sum()),
            "num_isolated_background": num_background - num_foreground,
            "num_union": num_background,
            "selected_id_order_sha256": ids_sha256(snps, architecture.background_indices),
            "all_foreground_nested_in_background": True,
            "fixed_across_replicates": True,
        },
        "h2_background_by_component": dict(zip(component_names, h2_background.tolist())),
        "omega_variance_by_component": dict(zip(component_names, omega_variances.tolist())),
        "expected_background_variance_by_component": dict(zip(component_names, h2_background.tolist())),
        "foreground_gamma_variance_ratio_to_background_variance": args.foreground_gamma_variance_ratio,
        "foreground_gamma_variance_by_component": {
            "GAMMA1": args.foreground_gamma_variance_ratio * omega_variances[1],
            "GAMMA2": args.foreground_gamma_variance_ratio * omega_variances[2],
        },
        "foreground_alpha_variance_ratio_to_background_variance": args.foreground_alpha_variance_ratio,
        "foreground_alpha_variance_by_component": {
            "ALPHA1": args.foreground_alpha_variance_ratio * omega_variances[0],
            "ALPHA2": args.foreground_alpha_variance_ratio * omega_variances[3],
        },
    }
    if args.paired_effect_streams:
        architecture_metadata["effect_stream_mode"] = "paired_by_original_snp"
        architecture_metadata["rng_streams"]["effects"] = {
            "background": "default_rng(SeedSequence([base_seed, 2, replicate, 0]))",
            "foreground": "default_rng(SeedSequence([base_seed, 2, replicate, 1]))",
            "indexing": "draw num_snps x 4 in original Step 2 SNP order, then subset",
        }
    save_checkpoint(architecture_frame, architecture_metadata, architecture_path, args.overwrite)
    atomic_write_table(
        pd.DataFrame(counts), out_dir / "architectures" / f"{scenario}_seed{args.seed}_counts.tsv.gz",
    )
    print(
        f"Fixed architecture: cores={len(core_indices)}, expanded foreground={num_foreground:,}, "
        f"background-only={num_background - num_foreground:,}; M_bg={num_background:,}; "
        f"omega variances={omega_variances.tolist()}", flush=True,
    )
    for chromosome, values in architecture.background_by_chromosome.items():
        print(
            f"  chr{chromosome}: foreground={values['num_foreground']:,}, "
            f"background-only={values['num_nonforeground_background']:,}, "
            f"buffer excluded={values['num_excluded_buffer_snps']:,}, "
            f"background total={values['num_background']:,}/{values['num_snps']:,}", flush=True,
        )

    architecture_sha256 = file_sha256(architecture_path)
    last_replicate = args.replicate_start + args.replicates - 1
    for replicate in range(args.replicate_start, last_replicate + 1):
        # Replicate-indexed streams make restarts and extensions order-independent.
        if args.paired_effect_streams:
            background_latent, foreground_latent = paired_replicate_latents(
                architecture, args.seed, replicate, len(snps),
            )
        else:
            effect_rng = np.random.default_rng(effect_seed(args.seed, replicate))
            background_latent, foreground_latent = replicate_latents(architecture, effect_rng)
        print(
            f"[{replicate}/{last_replicate}] Drawing effects: fixed foreground={num_foreground:,}, "
            f"background-only={num_background - num_foreground:,}; M_bg={num_background:,}", flush=True,
        )

        # 3. Reuse this replicate's latent draws across the setting grid.
        for rho in args.rhos:
            omega = background_covariance(rho, omega_variances)
            sigma_f = foreground_covariance(
                rho, omega_variances, args.foreground_gamma_variance_ratio, args.foreground_alpha_variance_ratio,
            )
            background = generate_background(architecture, background_latent, rho, omega_variances, len(snps))
            foreground = generate_foreground(
                architecture, foreground_latent, rho, omega_variances,
                args.foreground_gamma_variance_ratio, args.foreground_alpha_variance_ratio, len(snps),
            )
            qc = empirical_qc(snps, architecture, background, foreground)
            for beta in args.betas:
                output = build_output(snps, architecture, background, foreground, beta)
                for n_aux in args.n_aux:
                    for n_target in args.n_target:
                        output_path = output_path_for(
                            out_dir, target, args.chromosomes, rho, beta, n_aux, n_target,
                            args.foreground_gamma_variance_ratio, scenario, replicate,
                        )
                        metadata = {
                            **architecture_metadata,
                            "model": "XMR four-trait model with randomly selected fixed SNP architecture",
                            "replicate": replicate,
                            "source": XMR_SOURCE,
                            "architecture_file": str(architecture_path),
                            "architecture_file_sha256": architecture_sha256,
                            "populations": {"pop1_auxiliary": "WB", "pop2_target": target},
                            "setting_id": setting_id(beta, rho, n_aux, n_target),
                            "sample_sizes": {
                                "auxiliary_N1": n_aux, "target_N2": n_target,
                                "auxiliary_SE": math.sqrt(1 / n_aux), "target_SE": math.sqrt(1 / n_target),
                                "note": "Sample sizes are applied in Step 5; true effects do not depend on N",
                            },
                            "beta": beta, "beta1_auxiliary": beta, "beta2_target": beta,
                            "rho_background_u1_u2": rho, "rho_gamma1_gamma2": rho,
                            "Omega_order": list(component_names),
                            "R_Omega": background_correlation(rho).tolist(), "Omega": omega.tolist(),
                            "Sigma_f_order": ["ALPHA1", "GAMMA1", "GAMMA2", "ALPHA2"],
                            "Sigma_f": sigma_f.tolist(),
                            "C_order": ["Y1", "X1", "X2", "Y2"], "C": C_MATRIX.tolist(),
                            "true_effect_equations": {
                                "TRUE_X1": "U1 + GAMMA1", "TRUE_X2": "U2 + GAMMA2",
                                "TRUE_Y1": "V1 + ALPHA1 + beta1 * GAMMA1",
                                "TRUE_Y2": "V2 + ALPHA2 + beta2 * GAMMA2",
                            },
                            "contains_ld_projection": False, "contains_sampling_error": False, "qc": qc,
                        }
                        save_checkpoint(output, metadata, output_path, args.overwrite)
        print(f"[{replicate}/{last_replicate}] True effects completed for all requested settings.", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
