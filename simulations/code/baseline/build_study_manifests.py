#!/usr/bin/env python3
"""Build compact path manifests from completed simulation checkpoints."""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path


def float_label(value: float) -> str:
    return format(value, ".12g").replace("-", "m").replace(".", "p").replace("+", "")


def setting_id(beta: float, rho: float, n_aux: int, n_target: int) -> str:
    return (
        f"beta{float_label(beta)}_rho{float_label(rho)}_"
        f"naux{n_aux}_ntarget{n_target}"
    )


def atomic_tsv(rows: list[dict[str, object]], columns: list[str], path: Path) -> None:
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def build_truth_manifest(root: Path) -> None:
    rows: list[dict[str, object]] = []
    for metadata_path in sorted(
        (root / "true_effects").glob("*/xmr_truth_*.json")
    ):
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        data_path = metadata_path.with_suffix(".tsv.gz")
        if not data_path.is_file():
            raise FileNotFoundError(f"Missing truth table for {metadata_path}")
        beta = float(metadata["beta1_auxiliary"])
        rho = float(metadata["rho_background_u1_u2"])
        sample_sizes = metadata["sample_sizes"]
        n_aux = int(sample_sizes["auxiliary_N1"])
        n_target = int(sample_sizes["target_N2"])
        current_setting = setting_id(beta, rho, n_aux, n_target)
        if metadata_path.parent.name != current_setting:
            raise ValueError(f"Setting directory mismatch: {metadata_path}")
        rows.append(
            {
                "setting_id": current_setting,
                "replicate": int(metadata["replicate"]),
                "beta_true": beta,
                "rho": rho,
                "n_aux": n_aux,
                "n_target": n_target,
                "truth_file": str(data_path),
                "metadata_file": str(metadata_path),
            }
        )
    rows.sort(key=lambda row: (row["setting_id"], row["replicate"]))
    atomic_tsv(
        rows,
        [
            "setting_id",
            "replicate",
            "beta_true",
            "rho",
            "n_aux",
            "n_target",
            "truth_file",
            "metadata_file",
        ],
        root / "manifests/truth_manifest.tsv",
    )


def build_sumstats_manifest(root: Path) -> None:
    rows: list[dict[str, object]] = []
    for metadata_path in sorted(
        (root / "summary_statistics").glob("*/xmr_sumstats_*.json")
    ):
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        data_path = metadata_path.with_suffix(".tsv.gz")
        if not data_path.is_file():
            raise FileNotFoundError(f"Missing summary table for {metadata_path}")
        beta = float(metadata["beta"])
        rho = float(metadata["rho"])
        n_aux = int(metadata["auxiliary_sample_size"])
        n_target = int(metadata["sample_sizes"]["target_EAS"])
        current_setting = setting_id(beta, rho, n_aux, n_target)
        if metadata_path.parent.name != current_setting:
            raise ValueError(f"Setting directory mismatch: {metadata_path}")
        rows.append(
            {
                "setting_id": current_setting,
                "replicate": int(metadata["replicate"]),
                "beta_true": beta,
                "rho": rho,
                "n_aux": n_aux,
                "n_target": n_target,
                "sumstats_file": str(data_path),
                "metadata_file": str(metadata_path),
                "truth_file": metadata["input_truth_file"],
            }
        )
    rows.sort(key=lambda row: (row["setting_id"], row["replicate"]))
    atomic_tsv(
        rows,
        [
            "setting_id",
            "replicate",
            "beta_true",
            "rho",
            "n_aux",
            "n_target",
            "sumstats_file",
            "metadata_file",
            "truth_file",
        ],
        root / "manifests/sumstats_manifest.tsv",
    )


def build_background_manifest(root: Path) -> None:
    rows: list[dict[str, object]] = []
    pattern = re.compile(r"rep([0-9]+)$")
    for path in sorted((root / "xmr_background").glob("*/rep*/background_parameters.rds")):
        match = pattern.match(path.parent.name)
        if match is None:
            continue
        rows.append(
            {
                "setting_id": path.parent.parent.name,
                "replicate": int(match.group(1)),
                "background_file": str(path),
            }
        )
    rows.sort(key=lambda row: (row["setting_id"], row["replicate"]))
    atomic_tsv(
        rows,
        ["setting_id", "replicate", "background_file"],
        root / "manifests/xmr_background_manifest.tsv",
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.expanduser().resolve()
    (root / "manifests").mkdir(parents=True, exist_ok=True)
    build_truth_manifest(root)
    build_sumstats_manifest(root)
    build_background_manifest(root)
    print(f"Updated manifests below {root / 'manifests'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
