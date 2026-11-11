#!/usr/bin/env python3
"""Generate LD-aware four-trait XMR GWAS summary statistics.

For trait t in population p, this script simulates

    beta_hat_t = R_p @ beta_true_t + se_t * R_p^(1/2) @ z_t,

where each SNP's innovation vector is independently distributed as N(0, C).
White British/EUR LD is used for Y1/X1 and EAS LD for X2/Y2. Both chromosomes
2 and 19 are projected block by block. The reported standard errors remain
1/sqrt(N); diagonal inflation in C is therefore intentionally represented in
the LDSC intercept rather than absorbed into the reported standard errors.
"""

from __future__ import annotations

import os
import argparse
import glob
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
    normalized_chromosome,
    normalized_chromosomes,
    prepare_schema,
)

from step4_generate_xmr_true_effects import (
    ARCHITECTURE_VERSION, CORE_POPULATIONS, core_candidates,
    background_covariance, foreground_covariance,
)


DATA_ROOT = Path("simulations/data")
ROOT = DATA_ROOT / "reduced_confounding"
DEFAULT_LD_DIR = DATA_ROOT / "vilma_ld"
DEFAULT_OUT_DIR = ROOT / "summary_statistics"
DEFAULT_N_AUXILIARY_VALUES = (40_000, 60_000)
DEFAULT_N_TARGET = 10_000
DEFAULT_EUR_Y_H2 = 0.3
DEFAULT_EUR_X_H2 = 0.2
DEFAULT_EAS_X_H2 = 0.2
DEFAULT_EAS_Y_H2 = 0.3
DEFAULT_BACKGROUND_ONLY_FRACTION = 0.70
DEFAULT_N_TARGET_ONLY = 5
DEFAULT_N_AUX_ONLY = 5
DEFAULT_N_SHARED = 50
DEFAULT_N_ANCHORS = DEFAULT_N_SHARED
DEFAULT_GAMMA_VARIANCE_RATIO = 10.0
DEFAULT_ALPHA_VARIANCE_RATIO = 1.0

# Trait and C order: Y1, X1, X2, Y2.
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

SAMPLING_ERROR_SEED_RULE = "SeedSequence([base_seed, 3, replicate])"
TRUE_COLUMNS = ["TRUE_Y1", "TRUE_X1", "TRUE_X2", "TRUE_Y2"]
COMPONENT_COLUMNS = [
    "V1",
    "U1",
    "U2",
    "V2",
    "ALPHA1",
    "GAMMA1",
    "GAMMA2",
    "ALPHA2",
]
ANNOTATION_COLUMNS = [
    "XMR_ORDER",
    "EFFECT_CLASS",
    "FOREGROUND_SOURCE",
    "HAS_BACKGROUND",
    "Z_FOREGROUND",
    "IS_CORE_FOREGROUND",
    "LD_NEIGHBOR_OF_CORE",
    "IN_ARCHITECTURE_UNION",
    "IV_AUX",
    "IV_TARGET",
]
TRAIT_SUFFIXES = ["out.pop1", "exp.pop1", "exp.pop2", "out.pop2"]
TRAIT_LABELS = ["Y1", "X1", "X2", "Y2"]
TRAIT_POPULATIONS = ["EUR/White British", "EUR/White British", "EAS", "EAS"]
EUR_TRAITS = (0, 1)
EAS_TRAITS = (2, 3)


@dataclass(frozen=True)
class TruthInput:
    path: Path
    replicate: int
    metadata_path: Path
    beta: float
    rho: float
    n_aux: int
    n_target: int
    setting_id: str


@dataclass(frozen=True)
class SpectralBlock:
    population: str
    block_index: int
    global_indices: np.ndarray
    eigenvectors: np.ndarray
    eigenvalues: np.ndarray
    full_dimension: int
    minimum_raw_eigenvalue: float


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
    values: list[int] = []
    for value in parse_float_list(raw, option):
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
            "Generate chr2+chr19 LD-aware XMR summary statistics from Step 4 "
            "true effects."
        )
    )
    parser.add_argument("--ld-dir", type=Path, default=DEFAULT_LD_DIR)
    parser.add_argument(
        "--chromosomes",
        nargs="+",
        default=["2", "19"],
        help="Chromosomes separated by spaces or commas (default: 2 19)",
    )
    truth_group = parser.add_mutually_exclusive_group()
    truth_group.add_argument(
        "--truth-glob",
        default=None,
        help=(
            "Glob selecting Step 4 truth files; default reads the "
            "per-setting directories below <study-root>/true_effects"
        ),
    )
    truth_group.add_argument(
        "--truth-files",
        type=Path,
        nargs="+",
        help="Explicit Step 4 true-effect files",
    )
    parser.add_argument("--eur-schema", type=Path, default=None)
    parser.add_argument("--eas-schema", type=Path, default=None)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument(
        "--n-aux",
        nargs="+",
        default=[str(value) for value in DEFAULT_N_AUXILIARY_VALUES],
        help=(
            "One or more auxiliary sample sizes separated by spaces or commas "
            "(default: 40000 60000)"
        ),
    )
    parser.add_argument(
        "--betas",
        nargs="+",
        default=["0", "0.2"],
        help="Step 4 beta values (default: 0 0.2)",
    )
    parser.add_argument(
        "--rhos",
        nargs="+",
        default=["0.4", "0.8"],
        help="Step 4 rho values (default: 0.4 0.8)",
    )
    parser.add_argument(
        "--n-target",
        nargs="+",
        default=[str(DEFAULT_N_TARGET)],
        help="One or more target sample sizes separated by spaces or commas",
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
        "--h2-eur-y", type=float, default=None
    )
    parser.add_argument(
        "--h2-eur-x", type=float, default=None
    )
    parser.add_argument(
        "--h2-eas-x", type=float, default=None
    )
    parser.add_argument(
        "--h2-eas-y", type=float, default=None
    )
    parser.add_argument(
        "--background-only-fraction",
        type=float,
        default=DEFAULT_BACKGROUND_ONLY_FRACTION,
        help="Background-only fraction per chromosome; foreground has additional background effects",
    )
    parser.add_argument("--n-anchors", type=int, default=None,
                        help="Legacy shared-core count; defaults to --n-shared and must equal it")
    parser.add_argument(
        "--n-target-only", type=int, default=DEFAULT_N_TARGET_ONLY
    )
    parser.add_argument("--n-aux-only", type=int, default=DEFAULT_N_AUX_ONLY)
    parser.add_argument("--n-shared", type=int, default=DEFAULT_N_SHARED)
    parser.add_argument(
        "--foreground-gamma-variance-ratio",
        type=float,
        default=DEFAULT_GAMMA_VARIANCE_RATIO,
    )
    parser.add_argument(
        "--foreground-alpha-variance-ratio",
        type=float,
        default=DEFAULT_ALPHA_VARIANCE_RATIO,
    )
    parser.add_argument("--ldscore-quantiles", type=float, nargs=2, default=[0.45, 0.75])
    parser.add_argument("--core-kb", type=int, default=10000)
    parser.add_argument("--core-r2", type=float, default=0.001)
    parser.add_argument("--ld-kb", type=int, default=1000)
    parser.add_argument("--ld-r2", type=float, default=0.4)
    parser.add_argument("--background-exclusion-kb", type=int, default=1000)
    parser.add_argument("--background-exclusion-r2", type=float, default=0.001)
    parser.add_argument(
        "--expected-total-snps",
        type=int,
        default=0,
        help="Fail unless the combined Step 2 files have this many SNPs",
    )
    parser.add_argument(
        "--expected-replicates",
        type=int,
        default=50,
        help=(
            "Expected truth replicates per full setting; use 0 to disable"
        ),
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=5,
        help="Number of replicates projected through LD together",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=414,
        help="Master seed; replicate-specific GWAS-error streams are derived from it",
    )
    parser.add_argument(
        "--eigenvalue-tolerance",
        type=float,
        default=1e-12,
        help="Relative eigenvalue threshold used for LD square roots",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=20,
        help="Report progress after this many LD blocks; use 0 to disable",
    )
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.n_anchors is None:
        args.n_anchors = args.n_shared
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
    if any(abs(rho) >= 1 for rho in args.rhos):
        parser.error("Every --rhos value must satisfy -1 < rho < 1")
    h2_options = {
        "--h2-eur-y": args.h2_eur_y,
        "--h2-eur-x": args.h2_eur_x,
        "--h2-eas-x": args.h2_eas_x,
        "--h2-eas-y": args.h2_eas_y,
    }
    for option, value in h2_options.items():
        if not 0 < value < 1:
            parser.error(f"{option} must be in (0, 1)")
    if not 0 < args.background_only_fraction <= 1:
        parser.error("--background-only-fraction must be in (0, 1]")
    count_values = (
        args.n_anchors,
        args.n_target_only,
        args.n_aux_only,
        args.n_shared,
    )
    if any(value < 0 for value in count_values):
        parser.error("Foreground counts cannot be negative")
    if args.n_anchors > args.n_shared:
        parser.error("--n-anchors cannot exceed --n-shared")
    if args.n_anchors != args.n_shared:
        parser.error(
            "All shared cores are anchors; "
            "--n-anchors must equal --n-shared"
        )
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
    if args.expected_total_snps < 0 or args.expected_replicates < 0:
        parser.error("Expected counts cannot be negative")
    if args.batch_size < 1:
        parser.error("--batch-size must be at least 1")
    if not 0 <= args.seed <= 2_147_483_646:
        parser.error("--seed must be between 0 and 2147483646")
    if not 0 <= args.eigenvalue_tolerance < 1:
        parser.error("--eigenvalue-tolerance must satisfy 0 <= value < 1")
    if args.progress_every < 0:
        parser.error("--progress-every cannot be negative")
    return args


def file_sha256(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def ids_sha256(ids: Iterable[object]) -> str:
    text = "\n".join(str(value) for value in ids)
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def chromosome_label(chromosomes: list[str]) -> str:
    return "_".join(f"chr{chromosome}" for chromosome in chromosomes)


def float_label(value: float) -> str:
    return (
        format(value, ".12g")
        .replace("-", "m")
        .replace(".", "p")
        .replace("+", "")
    )


def setting_id(beta: float, rho: float, n_aux: int, n_target: int) -> str:
    return (
        f"beta{float_label(beta)}_rho{float_label(rho)}_"
        f"naux{n_aux}_ntarget{n_target}"
    )


def default_truth_globs(
    ld_dir: Path, args: argparse.Namespace
) -> list[str]:
    del ld_dir
    h2_label = (
        f"h2eury{float_label(args.h2_eur_y)}_"
        f"eurx{float_label(args.h2_eur_x)}_"
        f"easx{float_label(args.h2_eas_x)}_"
        f"easy{float_label(args.h2_eas_y)}"
    )
    return [
        str(
            ROOT
            / "true_effects"
            / "*"
            / (
                f"xmr_truth_*_{h2_label}_*"
                f"core{args.n_shared}s{args.n_target_only}t{args.n_aux_only}a_*_fixedIV_rep*.tsv.gz"
            )
        )
    ]


def companion_json(path: Path) -> Path:
    text = str(path)
    if not text.endswith(".tsv.gz"):
        raise ValueError(f"Truth file must end in .tsv.gz: {path}")
    return Path(text[: -len(".tsv.gz")] + ".json")


def replicate_from_path(path: Path, metadata_path: Path) -> int:
    try:
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"Invalid JSON in {metadata_path}: {exc}") from exc
    if isinstance(metadata, dict) and "replicate" in metadata:
        return int(metadata["replicate"])
    match = re.search(r"(?:^|_)rep(?:licate)?[_-]?(\d+)(?:_|\.|$)", path.name)
    if match is None:
        raise ValueError(f"Cannot determine replicate number from {path}")
    return int(match.group(1))


def truth_sample_sizes(
    metadata: dict[str, object],
    metadata_path: Path,
    args: argparse.Namespace,
) -> tuple[int, int]:
    sample_sizes = metadata.get("sample_sizes")
    if isinstance(sample_sizes, dict):
        n_aux = sample_sizes.get("auxiliary_N1")
        n_target = sample_sizes.get("target_N2")
        if n_aux is not None and n_target is not None:
            return int(n_aux), int(n_target)
    match = re.fullmatch(
        r"beta.+_rho.+_naux([0-9]+)_ntarget([0-9]+)",
        metadata_path.parent.name,
    )
    if match is not None:
        return int(match.group(1)), int(match.group(2))
    if len(args.n_aux) == 1 and len(args.n_target) == 1:
        return args.n_aux[0], args.n_target[0]
    raise ValueError(
        f"Cannot determine N_aux/N_target for {metadata_path}; regenerate "
        "Step 4 outputs in per-setting directories"
    )


def resolve_truth_inputs(args: argparse.Namespace, ld_dir: Path) -> list[TruthInput]:
    if args.truth_files is not None:
        paths = [path.expanduser().resolve() for path in args.truth_files]
        selector = "explicit --truth-files"
    elif args.truth_glob is not None:
        paths = [Path(value).resolve() for value in glob.glob(args.truth_glob)]
        selector = args.truth_glob
        if not paths:
            requested = Path(args.truth_glob)
            nested_glob = str(requested.parent / "*" / requested.name)
            paths = [
                Path(value).resolve() for value in glob.glob(nested_glob)
            ]
            if paths:
                print(
                    "No flat Step 4 truth files matched; using per-setting "
                    f"directories: {nested_glob}"
                )
                selector = nested_glob
            else:
                selector = f"{args.truth_glob}; fallback: {nested_glob}"
    else:
        patterns = default_truth_globs(ld_dir, args)
        paths = [
            Path(value).resolve()
            for pattern in patterns
            for value in glob.glob(pattern)
        ]
        selector = ", ".join(patterns)
    if not paths:
        raise FileNotFoundError(f"No Step 4 truth files matched: {selector}")
    inputs: list[TruthInput] = []
    for path in sorted(set(paths)):
        if not path.is_file():
            raise FileNotFoundError(f"Truth file does not exist: {path}")
        metadata_path = companion_json(path)
        if not metadata_path.is_file():
            raise FileNotFoundError(f"Missing truth metadata: {metadata_path}")
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        beta = float(metadata["beta1_auxiliary"])
        rho = float(metadata["rho_background_u1_u2"])
        n_aux, n_target = truth_sample_sizes(metadata, metadata_path, args)
        if not any(math.isclose(beta, value) for value in args.betas):
            continue
        if not any(math.isclose(rho, value) for value in args.rhos):
            continue
        if n_aux not in args.n_aux or n_target not in args.n_target:
            continue
        current_setting = setting_id(beta, rho, n_aux, n_target)
        observed_setting = metadata.get("setting_id")
        if observed_setting != current_setting or path.parent.name != current_setting:
            raise ValueError(
                f"Setting mismatch in {metadata_path}: metadata="
                f"{observed_setting!r}, directory={path.parent.name!r}, "
                f"expected={current_setting!r}"
            )
        inputs.append(
            TruthInput(
                path=path,
                replicate=replicate_from_path(path, metadata_path),
                metadata_path=metadata_path,
                beta=beta,
                rho=rho,
                n_aux=n_aux,
                n_target=n_target,
                setting_id=current_setting,
            )
        )
    if not inputs:
        raise FileNotFoundError(
            "Truth files were found, but none matched the requested setting grid"
        )
    inputs.sort(key=lambda item: (item.setting_id, item.replicate))
    keys = [(item.setting_id, item.replicate) for item in inputs]
    if len(set(keys)) != len(keys):
        raise ValueError("Selected truth files contain duplicate setting/replicate keys")
    for rho in args.rhos:
        for beta in args.betas:
            for n_aux in args.n_aux:
                for n_target in args.n_target:
                    selected = [
                        item
                        for item in inputs
                        if math.isclose(item.rho, rho)
                        and math.isclose(item.beta, beta)
                        and item.n_aux == n_aux
                        and item.n_target == n_target
                    ]
                    if (
                        args.expected_replicates
                        and len(selected) != args.expected_replicates
                    ):
                        raise ValueError(
                            f"{setting_id(beta, rho, n_aux, n_target)}: "
                            f"expected {args.expected_replicates:,} truth "
                            f"files but selected {len(selected):,}"
                        )
    return inputs


def load_step2_reference(
    ld_dir: Path,
    chromosomes: list[str],
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
        "L2_EUR",
        "L2_EAS",
        "L12_EUR_EAS",
    }
    for chromosome in chromosomes:
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
            raise ValueError(f"{path} is missing columns: {missing}")
        observed = {
            normalized_chromosome(value)
            for value in frame["CHROM"].dropna().unique()
        }
        if observed != {chromosome}:
            raise ValueError(
                f"{path} contains chromosomes {sorted(observed)}; "
                f"expected {chromosome}"
            )
        frame = frame.copy()
        frame["ID"] = frame["ID"].astype(str)
        frame["CHROM"] = chromosome
        frame["BP"] = pd.to_numeric(frame["BP"], errors="coerce")
        frame["A1"] = frame["A1"].astype(str).str.upper()
        frame["A2"] = frame["A2"].astype(str).str.upper()
        if not np.isfinite(frame["BP"].to_numpy(dtype=np.float64)).all():
            raise ValueError(f"Non-finite BP values in {path}")
        frame["BP"] = frame["BP"].astype(np.int64)
        score_columns = [
            column
            for column in frame.columns
            if column.startswith("L2") or column.startswith("L12")
        ]
        for column in score_columns:
            values = pd.to_numeric(frame[column], errors="coerce").to_numpy(
                dtype=np.float64
            )
            if not np.isfinite(values).all():
                raise ValueError(f"Non-finite values in {path}:{column}")
            frame[column] = values
        frames.append(frame)
        inputs.append(
            {
                "chromosome": chromosome,
                "path": str(path),
                "sha256": file_sha256(path),
                "num_snps": int(len(frame)),
            }
        )
    reference = pd.concat(frames, ignore_index=True, sort=False)
    duplicated = reference["ID"].duplicated(keep=False)
    if duplicated.any():
        examples = ", ".join(reference.loc[duplicated, "ID"].head(5))
        raise ValueError(f"Duplicate Step 2 SNP IDs; examples: {examples}")
    if expected_total_snps and len(reference) != expected_total_snps:
        raise ValueError(
            f"Expected {expected_total_snps:,} SNPs but found {len(reference):,}"
        )
    return reference, inputs


def rounded_fraction_count(total: int, fraction: float) -> int:
    return int(math.floor(total * fraction + 0.5))


def expected_background_counts(
    reference: pd.DataFrame,
    chromosomes: list[str],
    fraction: float,
    metadata_path: Path,
) -> dict[str, int]:
    # Every SNP in the union, including LD neighbors, has background effects.
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    selection = metadata.get("background_selection", {})
    if (
        selection.get("fraction_numerator")
        != "background_only"
        or selection.get("ld_neighbors_have_background_effects") is not True
        or selection.get("ld_neighbors_have_foreground_effects") is not True
    ):
        raise ValueError(f"{metadata_path}: old background design; rerun Step 4")
    values = reference["CHROM"].astype(str).to_numpy()
    counts = {}
    for chromosome in chromosomes:
        size = int(np.count_nonzero(values == chromosome))
        quota = rounded_fraction_count(size, fraction)
        info = selection.get("by_chromosome", {}).get(chromosome, {})
        raw_foreground = info.get("num_foreground", -1)
        foreground_count = int(raw_foreground)
        if foreground_count != raw_foreground or not 0 <= foreground_count <= size - quota:
            raise ValueError(f"{metadata_path}: invalid chr{chromosome} foreground count")
        counts[chromosome] = quota + foreground_count
        if (info.get("num_background") != counts[chromosome]
                or info.get("num_union") != counts[chromosome]
                or info.get("num_nonforeground_background") != quota):
            raise ValueError(f"{metadata_path}: chr{chromosome} background-only quota/M_bg disagree")
    if sum(counts.values()) <= 0:
        raise ValueError(f"{metadata_path}: no SNPs with direct background effects")
    return counts


def metadata_number(data: dict[str, object], key: str, path: Path) -> float:
    if key not in data:
        raise ValueError(f"{path} is missing metadata field {key!r}")
    try:
        value = float(data[key])
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{path} has invalid numeric field {key!r}") from exc
    if not math.isfinite(value):
        raise ValueError(f"{path} has non-finite metadata field {key!r}")
    return value


def metadata_number_mapping(
    data: dict[str, object], key: str, path: Path
) -> dict[str, float]:
    observed = data.get(key)
    if not isinstance(observed, dict):
        raise ValueError(f"{path} is missing metadata mapping {key!r}")
    result: dict[str, float] = {}
    for name, raw_value in observed.items():
        try:
            value = float(raw_value)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"{path} has invalid {key}[{name!r}]"
            ) from exc
        if not math.isfinite(value):
            raise ValueError(f"{path} has non-finite {key}[{name!r}]")
        result[str(name)] = value
    return result


def assert_close(
    observed: float,
    expected: float,
    description: str,
    path: Path,
    tolerance: float = 1e-12,
) -> None:
    if not math.isclose(observed, expected, rel_tol=1e-10, abs_tol=tolerance):
        raise ValueError(
            f"{path} has {description}={observed:.12g}; expected "
            f"{expected:.12g}"
        )


def load_and_validate_truth_metadata(
    truth: TruthInput,
    args: argparse.Namespace,
    reference: pd.DataFrame,
    step2_inputs: list[dict[str, object]],
    schema_paths: dict[str, Path],
    background_counts: dict[str, int],
) -> dict[str, object]:
    try:
        metadata = json.loads(truth.metadata_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"Invalid JSON in {truth.metadata_path}: {exc}") from exc
    if not isinstance(metadata, dict):
        raise ValueError(f"Truth metadata must be an object: {truth.metadata_path}")
    if int(metadata_number(metadata, "replicate", truth.metadata_path)) != truth.replicate:
        raise ValueError(f"Replicate mismatch in {truth.metadata_path}")
    if int(metadata_number(metadata, "num_snps", truth.metadata_path)) != len(reference):
        raise ValueError(f"SNP count mismatch in {truth.metadata_path}")
    if metadata.get("chromosomes") != args.chromosomes:
        raise ValueError(
            f"{truth.metadata_path} chromosomes differ from {args.chromosomes}"
        )
    if metadata.get("contains_ld_projection") is not False:
        raise ValueError(f"{truth.metadata_path} must contain pre-LD effects")
    if metadata.get("contains_sampling_error") is not False:
        raise ValueError(f"{truth.metadata_path} already contains sampling error")

    sample_sizes = metadata.get("sample_sizes")
    if not isinstance(sample_sizes, dict):
        raise ValueError(f"{truth.metadata_path} lacks sample_sizes metadata")
    if int(sample_sizes.get("auxiliary_N1", -1)) != truth.n_aux:
        raise ValueError(f"N_aux mismatch in {truth.metadata_path}")
    if int(sample_sizes.get("target_N2", -1)) != truth.n_target:
        raise ValueError(f"N_target mismatch in {truth.metadata_path}")
    if metadata.get("setting_id") != truth.setting_id:
        raise ValueError(f"Setting ID mismatch in {truth.metadata_path}")

    expected_h2 = {
        "V1": args.h2_eur_y,
        "U1": args.h2_eur_x,
        "U2": args.h2_eas_x,
        "V2": args.h2_eas_y,
    }
    observed_h2 = metadata_number_mapping(
        metadata, "h2_background_by_component", truth.metadata_path
    )
    if set(observed_h2) != set(expected_h2):
        raise ValueError(
            f"h2 background components differ in {truth.metadata_path}"
        )
    for component, expected in expected_h2.items():
        assert_close(
            observed_h2[component],
            expected,
            f"h2_background_by_component[{component}]",
            truth.metadata_path,
        )
    background_selection = metadata.get("background_selection")
    if not isinstance(background_selection, dict):
        raise ValueError(f"{truth.metadata_path} lacks background_selection")
    if metadata.get("architecture_version") != ARCHITECTURE_VERSION:
        raise ValueError(f"{truth.metadata_path}: regenerate Step 4 with the randomly selected fixed architecture")
    if metadata.get("paired_across_settings_within_replicate") is not True:
        raise ValueError(f"{truth.metadata_path}: settings must be paired within each replicate")
    if background_selection.get("fixed_across_replicates") is not True:
        raise ValueError(f"{truth.metadata_path}: background SNPs must be fixed across replicates")
    if background_selection.get("background_only_mutual_independence_required") is not False:
        raise ValueError(f"{truth.metadata_path}: background-only SNPs may have mutual LD")
    if background_selection.get("all_foreground_nested_in_background") is not True:
        raise ValueError(f"Foreground is not nested in background in {truth.metadata_path}")
    if (
        background_selection.get("fraction_numerator")
        != "background_only"
        or background_selection.get("M_bg_definition")
        != "expanded_foreground_plus_background_only"
        or background_selection.get("ld_neighbors_have_background_effects") is not True
        or background_selection.get("ld_neighbors_have_foreground_effects") is not True
        or background_selection.get("isolated_background_foreground_ld_independent") is not True
        or background_selection.get("background_isolation_reference") != "all_expanded_foreground"
    ):
        raise ValueError(f"{truth.metadata_path}: incompatible background/LD-neighbor design")
    neighborhood = background_selection.get("foreground_ld_neighborhood", {})
    if set(neighborhood.get("populations", [])) != {"EUR", "EAS"}:
        raise ValueError(f"{truth.metadata_path}: foreground LD neighborhoods must use EUR and EAS")
    assert_close(
        float(neighborhood.get("r2_threshold", math.nan)),
        args.ld_r2, "foreground LD-neighborhood r2", truth.metadata_path,
    )
    if neighborhood.get("window_kb") != args.ld_kb:
        raise ValueError(f"{truth.metadata_path}: foreground LD-neighborhood window mismatch")
    exclusion = background_selection.get("background_only_ld_exclusion", {})
    if (
        set(exclusion.get("populations", [])) != {"EUR", "EAS"}
        or exclusion.get("window_kb") != args.background_exclusion_kb
        or exclusion.get("r2_threshold") != args.background_exclusion_r2
    ):
        raise ValueError(f"{truth.metadata_path}: expanded-foreground exclusion settings mismatch")
    assert_close(
        float(background_selection.get("requested_background_only_fraction_per_chromosome", math.nan)),
        args.background_only_fraction,
        "background-only fraction",
        truth.metadata_path,
    )
    expected_total_background = sum(background_counts.values())
    if int(background_selection.get("num_background", -1)) != expected_total_background:
        raise ValueError(f"Background count mismatch in {truth.metadata_path}")
    by_chromosome = background_selection.get("by_chromosome")
    if not isinstance(by_chromosome, dict):
        raise ValueError(f"Missing chromosome background counts in {truth.metadata_path}")
    for chromosome, expected in background_counts.items():
        observed = by_chromosome.get(chromosome)
        if not isinstance(observed, dict) or int(observed.get("num_background", -1)) != expected:
            raise ValueError(
                f"Chromosome {chromosome} background count mismatch in "
                f"{truth.metadata_path}"
            )
    expected_omega = {
        component: value / expected_total_background
        for component, value in expected_h2.items()
    }
    observed_omega = metadata_number_mapping(
        metadata, "omega_variance_by_component", truth.metadata_path
    )
    if set(observed_omega) != set(expected_omega):
        raise ValueError(
            f"omega variance components differ in {truth.metadata_path}"
        )
    for component, expected in expected_omega.items():
        assert_close(
            observed_omega[component],
            expected,
            f"omega_variance_by_component[{component}]",
            truth.metadata_path,
            tolerance=1e-15,
        )
    assert_close(
        metadata_number(
            metadata,
            "foreground_gamma_variance_ratio_to_background_variance",
            truth.metadata_path,
        ),
        args.foreground_gamma_variance_ratio,
        "foreground gamma variance ratio",
        truth.metadata_path,
    )
    assert_close(
        metadata_number(
            metadata,
            "foreground_alpha_variance_ratio_to_background_variance",
            truth.metadata_path,
        ),
        args.foreground_alpha_variance_ratio,
        "foreground alpha variance ratio",
        truth.metadata_path,
    )
    expected_gamma_variance = {
        "GAMMA1": args.foreground_gamma_variance_ratio * expected_omega["U1"],
        "GAMMA2": args.foreground_gamma_variance_ratio * expected_omega["U2"],
    }
    observed_gamma_variance = metadata_number_mapping(
        metadata, "foreground_gamma_variance_by_component", truth.metadata_path
    )
    if set(observed_gamma_variance) != set(expected_gamma_variance):
        raise ValueError(
            f"foreground gamma components differ in {truth.metadata_path}"
        )
    for component, expected in expected_gamma_variance.items():
        assert_close(
            observed_gamma_variance[component],
            expected,
            f"foreground_gamma_variance_by_component[{component}]",
            truth.metadata_path,
            tolerance=1e-15,
        )
    expected_alpha_variance = {
        "ALPHA1": args.foreground_alpha_variance_ratio * expected_omega["V1"],
        "ALPHA2": args.foreground_alpha_variance_ratio * expected_omega["V2"],
    }
    observed_alpha_variance = metadata_number_mapping(
        metadata, "foreground_alpha_variance_by_component", truth.metadata_path
    )
    if set(observed_alpha_variance) != set(expected_alpha_variance):
        raise ValueError(
            f"foreground alpha components differ in {truth.metadata_path}"
        )
    for component, expected in expected_alpha_variance.items():
        assert_close(
            observed_alpha_variance[component],
            expected,
            f"foreground_alpha_variance_by_component[{component}]",
            truth.metadata_path,
            tolerance=1e-15,
        )

    foreground_selection = metadata.get("foreground_selection")
    class_assignment = metadata.get("foreground_class_assignment")
    if not isinstance(foreground_selection, dict) or not isinstance(
        class_assignment, dict
    ):
        raise ValueError(f"Missing foreground metadata in {truth.metadata_path}")
    n_core = args.n_target_only + args.n_aux_only + args.n_shared
    if foreground_selection.get("num_core") != n_core or foreground_selection.get("core_counts") != {
        "shared": args.n_shared, "target_only": args.n_target_only, "aux_only": args.n_aux_only,
    }:
        raise ValueError(f"{truth.metadata_path}: incorrect core counts")
    if foreground_selection.get("fixed_across_replicates") is not True:
        raise ValueError(f"{truth.metadata_path}: foreground SNPs must be fixed across replicates")
    if foreground_selection.get("selection_rule") != "random_order_greedy_LD_selection_without_replacement_no_retries":
        raise ValueError(f"{truth.metadata_path}: incorrect core selection rule")
    if (foreground_selection.get("independence_threshold_r2") != args.core_r2
            or foreground_selection.get("independence_window_kb") != args.core_kb):
        raise ValueError(f"{truth.metadata_path}: core independence settings mismatch")
    candidates, bounds = core_candidates(reference, args.ldscore_quantiles)
    if (foreground_selection.get("ldscore_quantiles") != args.ldscore_quantiles
            or foreground_selection.get("num_core_candidates") != len(candidates)
            or foreground_selection.get("ldscore_filter") != "intersection_of_both_inclusive_intervals_core_only"):
        raise ValueError(f"{truth.metadata_path}: core LD-score candidate settings mismatch")
    for population, expected in bounds.items():
        observed = foreground_selection.get("ldscore_bounds", {}).get(population, {})
        for bound in ("lower", "upper"):
            assert_close(float(observed.get(bound, math.nan)), expected[bound],
                         f"{population} LD-score {bound} quantile", truth.metadata_path)
    expansion = foreground_selection.get("ld_expansion", {})
    if (
        expansion.get("rule") != "one_hop_from_core_no_recursion"
        or expansion.get("neighbor_class_rule")
        != "union_of_linked_core_active_populations_core_classes_fixed"
        or expansion.get("effect_draws_across_snps")
        != "independent_same_per_snp_foreground_variance"
    ):
        raise ValueError(f"{truth.metadata_path}: incompatible foreground LD-expansion design")
    expanded_counts = expansion.get("expanded_class_counts", {})
    if set(expanded_counts) != {"shared", "target_only", "aux_only"}:
        raise ValueError(f"{truth.metadata_path}: invalid expanded foreground classes")
    for label, minimum in (("shared", args.n_shared), ("target_only", args.n_target_only),
                           ("aux_only", args.n_aux_only)):
        count = expanded_counts.get(label, -1)
        if not isinstance(count, int) or count < minimum or class_assignment.get(label) != count:
            raise ValueError(f"{truth.metadata_path}: invalid expanded {label} count")
    n_foreground = sum(expanded_counts.values())
    if (
        not n_core <= n_foreground <= expected_total_background
        or foreground_selection.get("num_foreground") != n_foreground
        or expansion.get("num_foreground") != n_foreground
        or expansion.get("num_core") != n_core
        or expansion.get("num_ld_neighbors") != n_foreground - n_core
    ):
        raise ValueError(f"{truth.metadata_path}: foreground core/neighbor totals disagree")
    if metadata.get("core_required_populations") != {
        label: list(populations) for label, populations in CORE_POPULATIONS.items()
    }:
        raise ValueError(f"{truth.metadata_path}: incorrect core population design")
    for label, expected in (("shared", args.n_shared), ("target_only", args.n_target_only),
                            ("aux_only", args.n_aux_only)):
        if class_assignment.get(f"core_{label}") != expected:
            raise ValueError(f"{truth.metadata_path}: incorrect {label} core count")
    if class_assignment.get("fixed_across_replicates") is not True:
        raise ValueError(f"{truth.metadata_path}: foreground classes must be fixed across replicates")
    if foreground_selection.get("groups_disjoint_by_snp_id") is not True:
        raise ValueError(f"Foreground SNP groups overlap in {truth.metadata_path}")

    metadata_c = np.asarray(metadata.get("C"), dtype=np.float64)
    if metadata_c.shape != C_MATRIX.shape or not np.allclose(
        metadata_c, C_MATRIX, rtol=0.0, atol=1e-12
    ):
        raise ValueError(f"C matrix mismatch in {truth.metadata_path}")
    beta1 = metadata_number(metadata, "beta1_auxiliary", truth.metadata_path)
    beta2 = metadata_number(metadata, "beta2_target", truth.metadata_path)
    assert_close(beta1, beta2, "beta1 versus beta2", truth.metadata_path)
    rho1 = metadata_number(metadata, "rho_background_u1_u2", truth.metadata_path)
    rho2 = metadata_number(metadata, "rho_gamma1_gamma2", truth.metadata_path)
    assert_close(rho1, rho2, "background versus foreground rho", truth.metadata_path)
    assert_close(beta1, truth.beta, "beta setting", truth.metadata_path)
    assert_close(rho1, truth.rho, "rho setting", truth.metadata_path)
    omega_variances = np.array([expected_omega[key] for key in ("V1", "U1", "U2", "V2")])
    for name, expected in (
        ("Omega", background_covariance(rho1, omega_variances)),
        ("Sigma_f", foreground_covariance(rho1, omega_variances,
                                          args.foreground_gamma_variance_ratio,
                                          args.foreground_alpha_variance_ratio)),
    ):
        observed = np.asarray(metadata.get(name), dtype=float)
        if observed.shape != expected.shape or not np.allclose(observed, expected, rtol=1e-10, atol=1e-15):
            raise ValueError(f"{truth.metadata_path}: incorrect {name}")

    metadata_step2 = metadata.get("step2_inputs")
    if not isinstance(metadata_step2, list):
        raise ValueError(f"Missing Step 2 input metadata in {truth.metadata_path}")
    observed_step2 = {
        str(item.get("chromosome")): item.get("sha256")
        for item in metadata_step2
        if isinstance(item, dict)
    }
    expected_step2 = {
        str(item["chromosome"]): item["sha256"] for item in step2_inputs
    }
    if observed_step2 != expected_step2:
        raise ValueError(f"Step 2 file hashes differ in {truth.metadata_path}")
    metadata_schemas = metadata.get("ld_schemas")
    if not isinstance(metadata_schemas, dict):
        raise ValueError(f"Missing LD schema metadata in {truth.metadata_path}")
    for population, schema_path in schema_paths.items():
        observed = metadata_schemas.get(population)
        if not isinstance(observed, dict) or observed.get("sha256") != file_sha256(
            schema_path
        ):
            raise ValueError(
                f"{population} schema hash differs in {truth.metadata_path}"
            )
    return metadata


def integer_column(frame: pd.DataFrame, column: str, path: Path) -> np.ndarray:
    values = pd.to_numeric(frame[column], errors="coerce").to_numpy(
        dtype=np.float64
    )
    if not np.isfinite(values).all() or not np.equal(values, np.rint(values)).all():
        raise ValueError(f"{path}:{column} must contain finite integers")
    return values.astype(np.int64)


def binary_column(frame: pd.DataFrame, column: str, path: Path) -> np.ndarray:
    values = integer_column(frame, column, path)
    if not np.isin(values, (0, 1)).all():
        raise ValueError(f"{path}:{column} must contain only 0/1")
    return values.astype(np.int8)


def validate_truth_design(
    truth: TruthInput,
    aligned: pd.DataFrame,
    effects: np.ndarray,
    components: np.ndarray,
    metadata: dict[str, object],
    args: argparse.Namespace,
    background_counts: dict[str, int],
) -> dict[str, object]:
    effect_class = aligned["EFFECT_CLASS"].astype(str).to_numpy()
    allowed = {"null_effect", "background", "target_only", "aux_only", "shared"}
    unexpected = sorted(set(effect_class) - allowed)
    if unexpected:
        raise ValueError(f"{truth.path} has unexpected classes: {unexpected}")
    num_background = sum(background_counts.values())
    expanded_counts = metadata["foreground_selection"]["ld_expansion"]["expanded_class_counts"]
    n_foreground = sum(expanded_counts.values())
    expected_counts = {
        **expanded_counts,
        "background": num_background - n_foreground,
        "null_effect": len(aligned) - num_background,
    }
    observed_counts = {
        str(key): int(value)
        for key, value in pd.Series(effect_class).value_counts().items()
    }
    if observed_counts != {key: value for key, value in expected_counts.items() if value > 0}:
        raise ValueError(
            f"{truth.path} class counts are {observed_counts}; expected "
            f"{expected_counts}"
        )
    has_background = binary_column(aligned, "HAS_BACKGROUND", truth.path)
    is_foreground = binary_column(aligned, "Z_FOREGROUND", truth.path)
    is_core = binary_column(aligned, "IS_CORE_FOREGROUND", truth.path)
    ld_neighbors = binary_column(aligned, "LD_NEIGHBOR_OF_CORE", truth.path)
    union = binary_column(aligned, "IN_ARCHITECTURE_UNION", truth.path)
    if np.any(is_core & ld_neighbors) or not np.array_equal(is_foreground, is_core | ld_neighbors):
        raise ValueError(f"{truth.path}: foreground must be the disjoint union of cores and LD neighbors")
    if np.any((is_foreground == 1) & (has_background != 1)):
        raise ValueError(f"{truth.path}: cores AND LD neighbors must have background effects")
    if not np.array_equal(union, has_background):
        raise ValueError(f"{truth.path}: every SNP in the union must have background effects")
    iv_aux = binary_column(aligned, "IV_AUX", truth.path)
    iv_target = binary_column(aligned, "IV_TARGET", truth.path)
    expected_background = (effect_class != "null_effect").astype(np.int8)
    expected_foreground = np.isin(
        effect_class, ("target_only", "aux_only", "shared")
    ).astype(np.int8)
    expected_iv_aux = np.isin(effect_class, ("aux_only", "shared")).astype(np.int8)
    expected_iv_target = np.isin(
        effect_class, ("target_only", "shared")
    ).astype(np.int8)
    for name, observed, expected in (
        ("HAS_BACKGROUND", has_background, expected_background),
        ("Z_FOREGROUND", is_foreground, expected_foreground),
        ("IV_AUX", iv_aux, expected_iv_aux),
        ("IV_TARGET", iv_target, expected_iv_target),
    ):
        if not np.array_equal(observed, expected):
            raise ValueError(f"{truth.path}:{name} disagrees with EFFECT_CLASS")
    chromosome_values = aligned["CHROM"].astype(str).to_numpy()
    for chromosome, expected in background_counts.items():
        chromosome_mask = chromosome_values == chromosome
        quota = rounded_fraction_count(int(chromosome_mask.sum()), args.background_only_fraction)
        neighbor_count = int(ld_neighbors[chromosome_mask].sum())
        chromosome_metadata = metadata["background_selection"]["by_chromosome"][chromosome]
        if (
            int(union[chromosome_mask].sum()) != expected
            or chromosome_metadata.get("num_union") != expected
            or int(np.count_nonzero(chromosome_mask & (effect_class == "background"))) != quota
            or chromosome_metadata.get("num_ld_neighbors") != neighbor_count
        ):
            raise ValueError(f"{truth.path}: chromosome {chromosome} background-only quota/LD-neighbor count mismatch")
        observed = int(
            np.count_nonzero(
                (chromosome_values == chromosome) & has_background.astype(bool)
            )
        )
        if observed != expected:
            raise ValueError(
                f"{truth.path} chromosome {chromosome} has {observed} "
                f"background SNPs; expected {expected}"
            )
    xmr_order = integer_column(aligned, "XMR_ORDER", truth.path)
    if not np.array_equal(
        np.sort(xmr_order), np.arange(1, len(aligned) + 1, dtype=np.int64)
    ):
        raise ValueError(f"{truth.path}:XMR_ORDER is not a full permutation")

    source = aligned["FOREGROUND_SOURCE"].astype(str).to_numpy()
    shared = source == "random_shared_core"
    target_only = source == "random_target_only_core"
    aux_only = source == "random_aux_only_core"
    if not np.array_equal(is_core.astype(bool), shared | target_only | aux_only):
        raise ValueError(f"{truth.path}: core annotations disagree with random-core sources")
    if not np.array_equal(ld_neighbors.astype(bool), source == "core_ld_neighbor"):
        raise ValueError(f"{truth.path}: LD-neighbor annotations disagree with foreground sources")
    if np.any(source[is_foreground == 0] != "not_foreground"):
        raise ValueError(f"{truth.path}: nonforeground SNPs have foreground sources")
    if int(shared.sum()) != args.n_shared:
        raise ValueError(f"{truth.path} has the wrong shared core count")
    if int(target_only.sum()) != args.n_target_only:
        raise ValueError(f"{truth.path} has the wrong target-only core count")
    if int(aux_only.sum()) != args.n_aux_only:
        raise ValueError(f"{truth.path} has the wrong auxiliary-only core count")
    if not np.all(effect_class[shared] == "shared"):
        raise ValueError(f"{truth.path} shared cores have incorrect classes")
    if not np.all(effect_class[target_only] == "target_only"):
        raise ValueError(f"{truth.path} target-only cores have incorrect classes")
    if not np.all(effect_class[aux_only] == "aux_only"):
        raise ValueError(f"{truth.path} auxiliary-only cores have incorrect classes")
    if not np.all((iv_aux[shared] == 1) & (iv_target[shared] == 1)):
        raise ValueError(f"{truth.path} shared SNPs are not IVs in both populations")
    if not np.all(iv_target[target_only] == 1) or not np.all(iv_aux[aux_only] == 1):
        raise ValueError(f"{truth.path} population-specific IV annotations are wrong")

    # Only cores must lie in both LD-score quantile intervals; neighbors may lie outside.
    for population, bounds in metadata["foreground_selection"]["ldscore_bounds"].items():
        scores = aligned[f"L2_{population}"].to_numpy(dtype=float)
        eligible = (scores >= bounds["lower"]) & (scores <= bounds["upper"])
        if not eligible[is_core.astype(bool)].all():
            raise ValueError(f"{truth.path}: cores outside {population} LD-score bounds")

    if not np.isfinite(components).all() or not np.isfinite(effects).all():
        raise ValueError(f"{truth.path} contains non-finite effects")
    if np.any(components[has_background == 0, :4] != 0.0):
        raise ValueError(f"{truth.path} has background effects outside background SNPs")
    if np.any(components[is_foreground == 0, 4:] != 0.0):
        raise ValueError(f"{truth.path} has foreground effects outside foreground SNPs")
    if np.any(components[effect_class == "target_only", 5] != 0.0):
        raise ValueError(f"{truth.path} target_only SNPs have nonzero GAMMA1")
    if np.any(components[effect_class == "aux_only", 6] != 0.0):
        raise ValueError(f"{truth.path} aux_only SNPs have nonzero GAMMA2")
    beta1 = float(metadata["beta1_auxiliary"])
    beta2 = float(metadata["beta2_target"])
    v1, u1, u2, v2, alpha1, gamma1, gamma2, alpha2 = components.T
    reconstructed = np.column_stack(
        (
            v1 + alpha1 + beta1 * gamma1,
            u1 + gamma1,
            u2 + gamma2,
            v2 + alpha2 + beta2 * gamma2,
        )
    )
    error = float(np.max(np.abs(effects - reconstructed)))
    if not np.allclose(effects, reconstructed, rtol=1e-10, atol=1e-12):
        raise ValueError(
            f"{truth.path} TRUE_* equations fail; max error={error:.6g}"
        )
    return {
        "effect_class_counts": observed_counts,
        "background_snps": int(has_background.sum()),
        "ld_neighbor_snps_with_background_effects": int(ld_neighbors.sum()),
        "ld_neighbor_snps_with_foreground_effects": int(ld_neighbors.sum()),
        "core_foreground_snps": int(is_core.sum()),
        "isolated_background_snps": int(
            np.count_nonzero((has_background == 1) & (is_foreground == 0) & (ld_neighbors == 0))
        ),
        "architecture_union_snps": int(union.sum()),
        "foreground_snps": int(is_foreground.sum()),
        "joint_shared_snps": int(shared.sum()),
        "eas_only_target_snps": int(target_only.sum()),
        "eur_only_auxiliary_snps": int(aux_only.sum()),
        "auxiliary_instruments": int(iv_aux.sum()),
        "target_instruments": int(iv_target.sum()),
        "selected_class_groups_disjoint_by_snp_id": bool(
            aligned.loc[shared | target_only | aux_only, "ID"].is_unique
        ),
        "all_foreground_nested_in_background": bool(
            np.all(has_background[is_foreground == 1] == 1)
        ),
        "shared_core_id_sha256": ids_sha256(aligned.loc[shared, "ID"]),
        "true_effect_equation_max_abs_error": error,
    }


def load_truth(
    truth: TruthInput,
    reference: pd.DataFrame,
    args: argparse.Namespace,
    step2_inputs: list[dict[str, object]],
    schema_paths: dict[str, Path],
    background_counts: dict[str, int],
) -> tuple[np.ndarray, pd.DataFrame, dict[str, object], dict[str, object]]:
    metadata = load_and_validate_truth_metadata(
        truth,
        args,
        reference,
        step2_inputs,
        schema_paths,
        background_counts,
    )
    header = pd.read_csv(truth.path, sep="\t", nrows=0)
    required = {
        "ID",
        "CHROM",
        "BP",
        "A1",
        "A2",
        *TRUE_COLUMNS,
        *COMPONENT_COLUMNS,
        *ANNOTATION_COLUMNS,
    }
    missing = sorted(required - set(header.columns))
    if missing:
        raise ValueError(f"{truth.path} is missing columns: {missing}")
    frame = pd.read_csv(
        truth.path,
        sep="\t",
        usecols=list(required),
        low_memory=False,
    )
    if len(frame) != len(reference):
        raise ValueError(f"{truth.path} has {len(frame)} rows; expected {len(reference)}")
    frame["ID"] = frame["ID"].astype(str)
    if frame["ID"].duplicated().any():
        raise ValueError(f"Duplicate IDs in {truth.path}")
    indexed = frame.set_index("ID")
    reference_ids = reference["ID"].astype(str)
    indexer = indexed.index.get_indexer(reference_ids)
    if (indexer < 0).any():
        examples = ", ".join(reference_ids[indexer < 0].head(5))
        raise ValueError(f"{truth.path} is missing SNPs; examples: {examples}")
    aligned = indexed.iloc[indexer].reset_index()
    observed_chromosome = aligned["CHROM"].map(normalized_chromosome).to_numpy()
    expected_chromosome = reference["CHROM"].astype(str).to_numpy()
    observed_bp = pd.to_numeric(aligned["BP"], errors="coerce").to_numpy(
        dtype=np.float64
    )
    expected_bp = reference["BP"].to_numpy(dtype=np.int64)
    observed_a1 = aligned["A1"].astype(str).str.upper().to_numpy()
    observed_a2 = aligned["A2"].astype(str).str.upper().to_numpy()
    mismatch = (
        (observed_chromosome != expected_chromosome)
        | (observed_bp != expected_bp)
        | (observed_a1 != reference["A1"].to_numpy())
        | (observed_a2 != reference["A2"].to_numpy())
    )
    if mismatch.any():
        raise ValueError(f"{truth.path} differs from Step 2 at {mismatch.sum()} SNPs")
    effects = aligned[TRUE_COLUMNS].apply(pd.to_numeric, errors="coerce").to_numpy(
        dtype=np.float64
    )
    components = aligned[COMPONENT_COLUMNS].apply(
        pd.to_numeric, errors="coerce"
    ).to_numpy(dtype=np.float64)
    aligned[["L2_EUR", "L2_EAS"]] = reference[["L2_EUR", "L2_EAS"]].to_numpy()
    annotations = aligned[ANNOTATION_COLUMNS].copy()
    design_qc = validate_truth_design(
        truth,
        aligned,
        effects,
        components,
        metadata,
        args,
        background_counts,
    )
    return effects, annotations, design_qc, metadata


def validate_fixed_assignment(
    annotations: pd.DataFrame, truth: TruthInput, reference_hash: str | None,
) -> str:
    digest = hashlib.sha256(
        pd.util.hash_pandas_object(annotations[ANNOTATION_COLUMNS], index=False).to_numpy().tobytes()
    ).hexdigest()
    if reference_hash is not None and digest != reference_hash:
        raise ValueError(f"Fixed SNP architecture differs across replicates/settings: {truth.path}")
    return digest


def spectral_decompose_schema(
    prepared: PreparedSchema,
    relative_tolerance: float,
    progress_every: int,
) -> tuple[list[SpectralBlock], dict[str, object]]:
    blocks: list[SpectralBlock] = []
    raw_minimum = math.inf
    total_dimension = 0
    retained_rank = 0
    clipped_negative = 0
    ordered = sorted(prepared.blocks.values(), key=lambda block: block.index)
    for count, source_block in enumerate(ordered, start=1):
        matrix = load_aligned_block(source_block)
        eigenvalues, eigenvectors = np.linalg.eigh(matrix)
        largest = float(max(np.max(np.abs(eigenvalues)), 1.0))
        tolerance = relative_tolerance * largest
        minimum = float(eigenvalues[0])
        raw_minimum = min(raw_minimum, minimum)
        if minimum < -max(tolerance, 1e-10):
            raise ValueError(
                f"{prepared.population} block {source_block.index} is not PSD: "
                f"minimum eigenvalue={minimum:.6g}"
            )
        clipped_negative += int(np.count_nonzero(eigenvalues < 0))
        keep = eigenvalues > tolerance
        if not keep.any():
            raise ValueError(
                f"{prepared.population} block {source_block.index} has zero rank"
            )
        values = np.asarray(eigenvalues[keep], dtype=np.float64)
        vectors = np.ascontiguousarray(eigenvectors[:, keep], dtype=np.float64)
        blocks.append(
            SpectralBlock(
                population=prepared.population,
                block_index=source_block.index,
                global_indices=source_block.global_indices.copy(),
                eigenvectors=vectors,
                eigenvalues=values,
                full_dimension=matrix.shape[0],
                minimum_raw_eigenvalue=minimum,
            )
        )
        total_dimension += matrix.shape[0]
        retained_rank += values.size
        if progress_every and (count % progress_every == 0 or count == len(ordered)):
            print(
                f"{prepared.population} spectral LD: processed "
                f"{count:,}/{len(ordered):,} blocks"
            )
    summary = {
        "population": prepared.population,
        "num_blocks": len(blocks),
        "total_dimension": total_dimension,
        "retained_rank": retained_rank,
        "retained_rank_fraction": retained_rank / total_dimension,
        "minimum_raw_eigenvalue": raw_minimum,
        "clipped_negative_eigenvalues": clipped_negative,
        "relative_eigenvalue_tolerance": relative_tolerance,
    }
    return blocks, summary


def project_blocks(
    source: np.ndarray,
    destination: np.ndarray,
    blocks: list[SpectralBlock],
    trait_indices: tuple[int, int],
    power: float,
) -> None:
    if source.shape != destination.shape or source.ndim != 3:
        raise ValueError("Source/destination must match [replicate,SNP,trait]")
    for block in blocks:
        indices = block.global_indices
        block_input = source[:, indices, :][:, :, trait_indices].copy()
        flat_input = block_input.transpose(1, 0, 2).reshape(len(indices), -1)
        coefficients = block.eigenvectors.T @ flat_input
        spectral_values = np.power(block.eigenvalues, power)
        projected = block.eigenvectors @ (spectral_values[:, None] * coefficients)
        projected = projected.reshape(
            len(indices), source.shape[0], len(trait_indices)
        ).transpose(1, 0, 2)
        for local_trait, global_trait in enumerate(trait_indices):
            destination[:, indices, global_trait] = projected[:, :, local_trait]


def validate_foreground_independence(
    annotations: pd.DataFrame,
    reference: pd.DataFrame,
    schemas: dict[str, PreparedSchema],
    args: argparse.Namespace,
) -> dict[str, object]:
    effect_class = annotations["EFFECT_CLASS"].astype(str).to_numpy()
    foreground_mask = np.isin(
        effect_class, ("shared", "target_only", "aux_only")
    )
    core_mask = annotations["IS_CORE_FOREGROUND"].to_numpy(dtype=bool)
    expected_count = args.n_target_only + args.n_aux_only + args.n_shared
    if int(core_mask.sum()) != expected_count:
        raise ValueError(f"Expected {expected_count} foreground core SNPs")
    bp = reference["BP"].to_numpy(dtype=np.int64)
    chromosome = reference["CHROM"].astype(str).to_numpy()
    sets = {
        "EUR_core_IVs": (
            np.flatnonzero(core_mask & np.isin(effect_class, ("shared", "aux_only"))),
            "EUR",
        ),
        "EAS_core_IVs": (
            np.flatnonzero(core_mask & np.isin(effect_class, ("shared", "target_only"))),
            "EAS",
        ),
    }
    set_results: dict[str, dict[str, int | float | bool]] = {}
    foreground_indices = np.flatnonzero(core_mask)
    expanded_indices = np.flatnonzero(foreground_mask)
    background_exclusion_mask = np.zeros(len(reference), dtype=bool)
    linked_by_class = {
        label: np.zeros(len(reference), dtype=bool)
        for label in ("shared", "target_only", "aux_only")
    }
    for set_name, (indices, population) in sets.items():
        prepared = schemas[population]
        pair_checks = 0
        pair_violations = 0
        maximum_pair_r2 = 0.0
        block_membership = prepared.membership[indices]
        # Expand only from cores, not recursively from their neighbors.
        all_foreground_membership = prepared.membership[foreground_indices]
        expanded_membership = prepared.membership[expanded_indices]
        for block_index in np.unique(expanded_membership):
            block = prepared.blocks[int(block_index)]
            block_foreground = indices[
                block_membership == block_index
            ]
            matrix = load_aligned_block(block)
            block_globals = block.global_indices
            all_block_foreground = foreground_indices[
                all_foreground_membership == block_index
            ]
            for index in all_block_foreground:
                in_window = (
                    (chromosome[block_globals] == chromosome[index])
                    & (np.abs(bp[block_globals] - bp[index]) <= args.ld_kb * 1000)
                )
                linked = in_window & (
                    np.square(matrix[block.local_by_global[int(index)], :])
                    >= args.ld_r2
                )
                linked_by_class[effect_class[index]][block_globals[linked]] = True
            for index in expanded_indices[expanded_membership == block_index]:
                in_window = (
                    (chromosome[block_globals] == chromosome[index])
                    & (np.abs(bp[block_globals] - bp[index]) <= args.background_exclusion_kb * 1000)
                )
                linked = in_window & (
                    np.square(matrix[block.local_by_global[int(index)], :])
                    >= args.background_exclusion_r2
                )
                background_exclusion_mask[block_globals[linked]] = True
            foreground_local = np.array(
                [block.local_by_global[int(index)] for index in block_foreground],
                dtype=np.int64,
            )
            for position, index in enumerate(block_foreground[:-1]):
                later = np.arange(position + 1, len(block_foreground))
                in_window = (
                    (chromosome[block_foreground[later]] == chromosome[index])
                    & (np.abs(bp[block_foreground[later]] - bp[index]) <= args.core_kb * 1000)
                )
                later = later[in_window]
                if not later.size:
                    continue
                values = np.square(
                    matrix[foreground_local[position], foreground_local[later]]
                )
                pair_checks += int(values.size)
                maximum_pair_r2 = max(maximum_pair_r2, float(np.max(values)))
                pair_violations += int(
                    np.count_nonzero(values >= args.core_r2)
                )
        if pair_violations:
            raise ValueError(
                f"{set_name} independence failed in {population}: "
                f"{pair_violations} pairs have r^2 >= "
                f"{args.core_r2:g} within "
                f"{args.core_kb / 1000:g} Mb"
            )
        set_results[set_name] = {
            "population": population,
            "num_snps": int(len(indices)),
            "pairs_checked_within_window": pair_checks,
            "pair_violations": pair_violations,
            "maximum_pair_r2": maximum_pair_r2,
            "passed": True,
        }
    for mask in linked_by_class.values():
        mask[foreground_indices] = False
    active_aux = linked_by_class["shared"] | linked_by_class["aux_only"]
    active_target = linked_by_class["shared"] | linked_by_class["target_only"]
    observed_ld_neighbors = active_aux | active_target
    annotated_neighbors = annotations["LD_NEIGHBOR_OF_CORE"].to_numpy(dtype=bool)
    if not np.array_equal(observed_ld_neighbors, annotated_neighbors):
        raise ValueError("LD-neighbor annotations disagree with the EUR/EAS LD matrices")
    if not np.array_equal(foreground_mask, core_mask | observed_ld_neighbors):
        raise ValueError("Foreground support differs from cores plus their one-hop LD neighbors")
    expected_labels = {
        "shared": active_aux & active_target,
        "target_only": active_target & ~active_aux,
        "aux_only": active_aux & ~active_target,
    }
    for label, selected in expected_labels.items():
        if not np.all(effect_class[selected] == label):
            raise ValueError(f"LD-neighbor {label} labels disagree with linked core classes")
    isolated_background = effect_class == "background"
    if np.any(background_exclusion_mask & isolated_background):
        raise ValueError("Background-only SNPs overlap the LD neighborhood of expanded foreground")
    return {
        "num_foreground": int(foreground_mask.sum()),
        "num_core": int(core_mask.sum()),
        "groups_disjoint_by_snp_id": True,
        "independence_rule": (
            "EUR shared+aux_only cores and EAS shared+target_only cores each "
            "have raw r^2 < threshold; expanded foreground neighbors may be correlated"
        ),
        "r2_threshold": args.core_r2,
        "window_kb": args.core_kb,
        "expansion_r2": args.ld_r2, "expansion_kb": args.ld_kb,
        "background_exclusion_r2": args.background_exclusion_r2,
        "background_exclusion_kb": args.background_exclusion_kb,
        "background_only_mutual_independence_required": False,
        "active_iv_sets": set_results,
        "isolated_background_independent_of_all_foreground_in_EUR_and_EAS": True,
        "num_excluded_buffer_snps": int(np.count_nonzero(background_exclusion_mask & ~foreground_mask)),
        "num_verified_ld_neighbors_with_background": int(observed_ld_neighbors.sum()),
        "num_verified_ld_neighbors_with_foreground": int(observed_ld_neighbors.sum()),
        "cross_block_correlations_assumed_zero": True,
    }


def noise_seed(base_seed: int, replicate: int) -> np.random.SeedSequence:
    return np.random.SeedSequence([base_seed, 3, replicate])


def project_batch(
    effects: np.ndarray,
    replicates: list[int],
    eur_blocks: list[SpectralBlock],
    eas_blocks: list[SpectralBlock],
    base_seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    num_replicates, num_snps, num_traits = effects.shape
    if num_traits != 4 or num_replicates != len(replicates):
        raise ValueError("Unexpected true-effect batch shape")
    marginal_means = np.zeros_like(effects)
    project_blocks(effects, marginal_means, eur_blocks, EUR_TRAITS, power=1.0)
    project_blocks(effects, marginal_means, eas_blocks, EAS_TRAITS, power=1.0)
    c_sqrt = np.linalg.cholesky(C_MATRIX)
    standardized_errors = np.empty_like(effects)
    for batch_index, replicate in enumerate(replicates):
        rng = np.random.default_rng(noise_seed(base_seed, replicate))
        standardized_errors[batch_index] = (
            rng.standard_normal((num_snps, 4)) @ c_sqrt.T
        )
    project_blocks(
        standardized_errors,
        standardized_errors,
        eur_blocks,
        EUR_TRAITS,
        power=0.5,
    )
    project_blocks(
        standardized_errors,
        standardized_errors,
        eas_blocks,
        EAS_TRAITS,
        power=0.5,
    )
    if not np.isfinite(marginal_means).all() or not np.isfinite(
        standardized_errors
    ).all():
        raise ValueError("LD projections contain non-finite values")
    return marginal_means, standardized_errors


def two_sided_normal_p(estimates: np.ndarray, standard_error: float) -> np.ndarray:
    scaled = np.abs(np.asarray(estimates, dtype=np.float64) / standard_error)
    return np.fromiter(
        (math.erfc(value / math.sqrt(2.0)) for value in scaled),
        dtype=np.float64,
        count=scaled.size,
    )


def numeric_summary(values: np.ndarray) -> dict[str, float]:
    return {
        "min": float(np.min(values)),
        "mean": float(np.mean(values)),
        "sd": float(np.std(values, ddof=1)),
        "max": float(np.max(values)),
    }


def build_output(
    reference: pd.DataFrame,
    annotations: pd.DataFrame,
    true_effects: np.ndarray,
    beta_hat: np.ndarray,
    sample_sizes: np.ndarray,
    standard_errors: np.ndarray,
) -> pd.DataFrame:
    if true_effects.shape != beta_hat.shape or true_effects.shape != (
        len(reference),
        4,
    ):
        raise ValueError("True effects and beta_hat must have shape [SNP,4]")
    output = pd.DataFrame(
        {
            "SNP": reference["ID"].astype(str),
            "ID": reference["ID"].astype(str),
            "CHR": reference["CHROM"],
            "BP": reference["BP"],
            "A1": reference["A1"],
            "A2": reference["A2"],
        }
    )
    for column in annotations.columns:
        output[column] = annotations[column].to_numpy()
    for trait_index, column in enumerate(TRUE_COLUMNS):
        output[column] = true_effects[:, trait_index]
    for trait_index, suffix in enumerate(TRAIT_SUFFIXES):
        standard_error = float(standard_errors[trait_index])
        estimates = beta_hat[:, trait_index]
        output[f"b.{suffix}"] = estimates
        output[f"se.{suffix}"] = standard_error
        output[f"pval.{suffix}"] = two_sided_normal_p(
            estimates, standard_error
        )
        output[f"N.{suffix}"] = int(sample_sizes[trait_index])
    score_columns = [
        column
        for column in reference.columns
        if column.startswith("L2") or column.startswith("L12")
    ]
    for column in score_columns:
        output[column] = reference[column].to_numpy()
    if "L2.pop1" not in output.columns:
        output["L2.pop1"] = reference["L2_EUR"].to_numpy()
    if "L2.pop2" not in output.columns:
        output["L2.pop2"] = reference["L2_EAS"].to_numpy()
    if "L12" not in output.columns:
        output["L12"] = reference["L12_EUR_EAS"].to_numpy()
    return output


def output_paths(truth: TruthInput, out_dir: Path) -> tuple[Path, Path]:
    stem = truth.path.name[: -len(".tsv.gz")]
    if stem.startswith("xmr_truth_"):
        stem = "xmr_sumstats_" + stem[len("xmr_truth_") :]
    else:
        stem += "_sumstats"
    if re.search(r"_rep([0-9]+)$", stem) is None:
        raise ValueError(f"Truth filename lacks a replicate suffix: {truth.path}")
    setting_dir = out_dir / truth.setting_id
    return setting_dir / f"{stem}.tsv.gz", setting_dir / f"{stem}.json"


def preflight_outputs(
    truths: list[TruthInput],
    out_dir: Path,
    overwrite: bool,
) -> None:
    if overwrite:
        return
    for truth in truths:
        data_path, metadata_path = output_paths(truth, out_dir)
        if data_path.exists() != metadata_path.exists():
            raise FileExistsError(
                "Incomplete checkpoint; remove it or use --overwrite:\n  "
                f"{data_path}\n  {metadata_path}"
            )


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


def atomic_write_table(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    frame.to_csv(
        temporary,
        sep="\t",
        index=False,
        float_format="%.12g",
        compression="gzip",
    )
    temporary.replace(path)


def atomic_write_json(data: dict[str, object], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary.replace(path)


def truth_scenario_summary(metadata: dict[str, object]) -> dict[str, object]:
    return {
        key: metadata[key]
        for key in (
            "architecture_version",
            "architecture_id",
            "architecture_file",
            "paired_across_settings_within_replicate",
            "chromosomes",
            "h2_background_by_component",
            "omega_variance_by_component",
            "foreground_gamma_variance_ratio_to_background_variance",
            "foreground_gamma_variance_by_component",
            "foreground_alpha_variance_ratio_to_background_variance",
            "foreground_alpha_variance_by_component",
            "sample_sizes",
            "foreground_selection",
            "foreground_class_assignment",
            "background_selection",
            "beta1_auxiliary",
            "beta2_target",
            "rho_background_u1_u2",
            "rho_gamma1_gamma2",
        )
    }


def main() -> int:
    args = parse_args()
    ld_dir = args.ld_dir.expanduser().resolve()
    if not ld_dir.is_dir():
        raise FileNotFoundError(f"LD directory does not exist: {ld_dir}")
    out_dir = args.out_dir.expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
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
    for population, path in schema_paths.items():
        if not path.is_file():
            raise FileNotFoundError(f"{population} schema does not exist: {path}")
    truths = resolve_truth_inputs(args, ld_dir)
    preflight_outputs(truths, out_dir, args.overwrite)
    reference, step2_inputs = load_step2_reference(
        ld_dir, args.chromosomes, args.expected_total_snps
    )
    eigenvalues_c = np.linalg.eigvalsh(C_MATRIX)
    if eigenvalues_c[0] <= 0:
        raise ValueError("C_MATRIX must be positive definite")

    prepared = {
        population: prepare_schema(path, reference, population)
        for population, path in schema_paths.items()
    }
    eur_blocks, eur_spectral_qc = spectral_decompose_schema(
        prepared["EUR"], args.eigenvalue_tolerance, args.progress_every
    )
    eas_blocks, eas_spectral_qc = spectral_decompose_schema(
        prepared["EAS"], args.eigenvalue_tolerance, args.progress_every
    )
    print(f"Selected {len(truths):,} Step 4 setting/replicate files")
    print(
        f"SNPs={len(reference):,}, background-only fraction={args.background_only_fraction:g}; "
        f"auxiliary N={','.join(str(value) for value in args.n_aux)}, "
        f"target N={','.join(str(value) for value in args.n_target)}"
    )

    shared_metadata = {
        "model": "LD-aware four-trait XMR summary-statistics simulation",
        "chromosomes": args.chromosomes,
        "trait_order": TRAIT_LABELS,
        "trait_populations": TRAIT_POPULATIONS,
        "C": C_MATRIX.tolist(),
        "C_order": TRAIT_LABELS,
        "C_eigenvalues": eigenvalues_c.tolist(),
        "formula": (
            "beta_hat_t = R_population @ beta_true_t + "
            "SE_t * R_population^(1/2) @ z_t"
        ),
        "latent_error": "for each SNP, z=(z_Y1,z_X1,z_X2,z_Y2) ~ N(0,C)",
        "reported_se_interpretation": (
            "SE_t=1/sqrt(N_t); C diagonal is retained as sampling-error/"
            "LDSC-intercept inflation and is not absorbed into reported SE"
        ),
        "sampling_error_covariance": (
            "within a population: SE_t*SE_s*C[t,s]*R_population; across "
            "populations: SE_t*SE_s*C[t,s]*R_EUR^(1/2)@R_EAS^(1/2).T"
        ),
        "background_h2_by_component": {
            "V1": args.h2_eur_y,
            "U1": args.h2_eur_x,
            "U2": args.h2_eas_x,
            "V2": args.h2_eas_y,
        },
        "background_only_fraction": args.background_only_fraction,
        "M_bg_definition": "expanded_foreground_plus_background_only",
        "foreground_design": {
            "core_target_only": args.n_target_only,
            "core_aux_only": args.n_aux_only,
            "core_shared": args.n_shared,
            "foreground_support": "core_plus_one_hop_EUR_EAS_LD_neighbors",
            "selection": "random_order_greedy_from_both_LDscore_quantile_intervals",
            "ldscore_quantiles": args.ldscore_quantiles,
            "fixed_across_replicates": True,
            "expansion_r2": args.ld_r2, "expansion_kb": args.ld_kb,
            "background_exclusion_r2": args.background_exclusion_r2,
            "background_exclusion_kb": args.background_exclusion_kb,
            "groups_disjoint_by_snp_id": True,
            "core_iv_independence_r2": args.core_r2,
            "core_iv_independence_kb": args.core_kb,
        },
        "step2_inputs": step2_inputs,
        "schemas": {
            population: {
                "path": str(path),
                "sha256": file_sha256(path),
            }
            for population, path in schema_paths.items()
        },
        "eur_spectral_qc": eur_spectral_qc,
        "eas_spectral_qc": eas_spectral_qc,
        "eigenvalue_tolerance": args.eigenvalue_tolerance,
        "master_seed": args.seed,
        "sampling_error_base_seed": args.seed,
        "sampling_error_seed_rule": SAMPLING_ERROR_SEED_RULE,
    }

    architecture_hash: str | None = None
    foreground_ld_qc: dict[str, object] | None = None
    for batch_start in range(0, len(truths), args.batch_size):
        batch_truths = truths[batch_start : batch_start + args.batch_size]
        effects_list: list[np.ndarray] = []
        annotations_list: list[pd.DataFrame] = []
        design_qc_list: list[dict[str, object]] = []
        metadata_list: list[dict[str, object]] = []
        for truth in batch_truths:
            background_counts = expected_background_counts(
                reference, args.chromosomes, args.background_only_fraction, truth.metadata_path,
            )
            effects, annotations, design_qc, metadata = load_truth(
                truth,
                reference,
                args,
                step2_inputs,
                schema_paths,
                background_counts,
            )
            architecture_hash = validate_fixed_assignment(annotations, truth, architecture_hash)
            if foreground_ld_qc is None:
                foreground_ld_qc = validate_foreground_independence(
                    annotations, reference, prepared, args,
                )
                print(
                    f"Verified fixed SNP architecture at replicate {truth.replicate:03d}: "
                    f"foreground={design_qc['foreground_snps']:,}; "
                    "core independence, one-hop expansion and background-only isolation passed.",
                    flush=True,
                )
            design_qc["foreground_ld_qc"] = foreground_ld_qc
            effects_list.append(effects)
            annotations_list.append(annotations)
            design_qc_list.append(design_qc)
            metadata_list.append(metadata)
        effect_batch = np.stack(effects_list, axis=0)
        replicates = [truth.replicate for truth in batch_truths]
        marginal_means, standardized_errors = project_batch(
            effect_batch,
            replicates,
            eur_blocks,
            eas_blocks,
            args.seed,
        )
        for batch_index, truth in enumerate(batch_truths):
            sample_sizes = np.array(
                [truth.n_aux, truth.n_aux, truth.n_target, truth.n_target],
                dtype=np.int64,
            )
            standard_errors = np.sqrt(
                1.0 / sample_sizes.astype(np.float64)
            )
            sampling_error = (
                standardized_errors[batch_index] * standard_errors[None, :]
            )
            beta_hat = marginal_means[batch_index] + sampling_error
            if not np.isfinite(beta_hat).all():
                raise ValueError(
                    "Generated summary statistics contain non-finite values"
                )
            sample_metadata = {
                "master_seed": args.seed,
                "sample_sizes": {
                    "auxiliary_EUR_WB": truth.n_aux,
                    "target_EAS": truth.n_target,
                },
                "reported_standard_errors": {
                    "auxiliary_EUR_WB": float(standard_errors[0]),
                    "target_EAS": float(standard_errors[2]),
                },
            }
            output_path, metadata_path = output_paths(truth, out_dir)
            if output_pair_is_complete(
                output_path, metadata_path, args.overwrite
            ):
                saved = json.loads(metadata_path.read_text(encoding="utf-8"))
                expected_checkpoint = {
                    "input_truth_file_sha256": file_sha256(truth.path),
                    "input_truth_metadata_sha256": file_sha256(truth.metadata_path),
                    "sampling_error_base_seed": args.seed,
                    "sampling_error_seed_rule": SAMPLING_ERROR_SEED_RULE,
                    "eigenvalue_tolerance": args.eigenvalue_tolerance,
                }
                for key, expected in expected_checkpoint.items():
                    if saved.get(key) != expected:
                        raise ValueError(f"{output_path}: changed {key}; use --overwrite or a new work root")
                print(
                    f"Skipped {truth.setting_id}, "
                    f"replicate {truth.replicate:03d}: "
                    f"{output_path}"
                )
                continue
            output = build_output(
                reference,
                annotations_list[batch_index],
                effect_batch[batch_index],
                beta_hat,
                sample_sizes,
                standard_errors,
            )
            trait_qc: dict[str, object] = {}
            for trait_index, label in enumerate(TRAIT_LABELS):
                trait_qc[label] = {
                    "true_effect": numeric_summary(
                        effect_batch[batch_index, :, trait_index]
                    ),
                    "marginal_mean": numeric_summary(
                        marginal_means[batch_index, :, trait_index]
                    ),
                    "sampling_error": numeric_summary(
                        sampling_error[:, trait_index]
                    ),
                    "beta_hat": numeric_summary(beta_hat[:, trait_index]),
                }
            result_metadata = dict(shared_metadata)
            current_truth = metadata_list[batch_index]
            result_metadata.update({
                "background_snps": current_truth["background_selection"]["num_background"],
                "omega_variance_by_component": current_truth["omega_variance_by_component"],
                "background_counts_by_chromosome": {
                    chrom: values["num_background"]
                    for chrom, values in current_truth["background_selection"]["by_chromosome"].items()
                },
            })
            result_metadata.update(sample_metadata)
            result_metadata.update(
                {
                    "setting_id": truth.setting_id,
                    "replicate": truth.replicate,
                    "beta": truth.beta,
                    "rho": truth.rho,
                    "auxiliary_sample_size": truth.n_aux,
                    "target_sample_size": truth.n_target,
                    "input_truth_file": str(truth.path),
                    "input_truth_file_sha256": file_sha256(truth.path),
                    "input_truth_metadata": str(truth.metadata_path),
                    "input_truth_metadata_sha256": file_sha256(
                        truth.metadata_path
                    ),
                    "truth_scenario": truth_scenario_summary(
                        metadata_list[batch_index]
                    ),
                    "truth_design_qc": design_qc_list[batch_index],
                    "num_snps": len(reference),
                    "output_file": str(output_path),
                    "trait_qc": trait_qc,
                }
            )
            atomic_write_table(output, output_path)
            atomic_write_json(result_metadata, metadata_path)
            print(
                f"Wrote {truth.setting_id}, "
                f"replicate {truth.replicate:03d}: {output_path}"
            )
        del effect_batch, marginal_means, standardized_errors
        print(
            f"Completed {min(batch_start + len(batch_truths), len(truths)):,}/"
            f"{len(truths):,} replicates"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
