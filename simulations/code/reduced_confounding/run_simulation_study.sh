#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_ROOT="simulations/data"
# C order: EUR_Y, EUR_X, EAS_X, EAS_Y; both cross-population blocks are zero.
WORK_ROOT="${DATA_ROOT}/reduced_confounding"
STAGE="all"
METHOD="all"
BETAS=(0 0.2)
RHOS=(0.4 0.8)
N_AUX_VALUES=(40000 60000)
N_TARGET_VALUES=(10000)
DEFAULT_N_TARGET_ONLY=5
DEFAULT_N_AUX_ONLY=5
DEFAULT_N_SHARED=50
N_TARGET_ONLY="${N_TARGET_ONLY-${DEFAULT_N_TARGET_ONLY}}"
N_AUX_ONLY="${N_AUX_ONLY-${DEFAULT_N_AUX_ONLY}}"
N_SHARED="${N_SHARED-${DEFAULT_N_SHARED}}"  # legacy anchor arguments use this value
WORK_ROOT_EXPLICIT=0
BACKGROUND_ONLY_FRACTION=0.60  # foreground adds extra background-effect SNPs
EUR_Y_H2=0.2
EUR_X_H2=0.2
EAS_X_H2=0.2
EAS_Y_H2=0.2
FOREGROUND_GAMMA_VARIANCE_RATIO=10
FOREGROUND_ALPHA_VARIANCE_RATIO=1
CORE_LD_R2=0.001
CORE_LD_KB=10000
FOREGROUND_LD_R2=0.8
FOREGROUND_LD_KB=1000
DEFAULT_BACKGROUND_LD_R2=0.001
BACKGROUND_LD_R2="${BACKGROUND_LD_R2-${DEFAULT_BACKGROUND_LD_R2}}"
BACKGROUND_LD_KB=1000
LDSCORE_QUANTILES=(0.45 0.85)
REPLICATES=50
THRESHOLDS="5e-5,5e-8"
PYTHON="${SIM_PYTHON:-python3}"
OVERWRITE=0
NO_FEMA_QC=1
XMR_PRIORITY="FEMA_P"
MASTER_SEED=1234
PAIRED_EFFECT_STREAMS=0

usage() {
  cat <<'EOF'
Usage:
  bash run_simulation_study.sh [options]

Options:
  --stage step3|data|background|methods|evaluate|all|check
  --method XMR|CAMeRa|MRBEE|TEMR|MRAPSS|all
                       Select methods for all/methods/evaluate (default: all)
  --beta VALUE...       One value or a comma/space-separated set
  --rho VALUE...        One value or a comma/space-separated set
  --n-aux VALUE...      One value or a comma/space-separated set
  --n-target VALUE...   One value or a comma/space-separated set
  --n-shared N          Shared core SNP count before LD expansion (positive integer; default: 50)
  --n-target-only N     Target-only core SNP count (nonnegative integer; default: 5)
  --n-aux-only N        Auxiliary-only core SNP count (nonnegative integer; default: 5)
                       Env: N_SHARED, N_TARGET_ONLY, N_AUX_ONLY; CLI takes precedence
  --rhos VALUES        Foreground AND background EUR-X/EAS-X correlations (coupled)
  --background-eurx-easy-correlation VALUE  Background EUR X/EAS Y correlation (default: 0.3)
  --background-eury-easx-correlation VALUE  Background EUR Y/EAS X correlation (default: 0.1)
  --eas-background-correlation VALUE  EAS background X/Y correlation (default: 0.4)
  --eas-c11 VALUE       EAS X intercept (default: 1.10)
  --eas-c12 VALUE       EAS X/Y intercept (default: 0.15)
  --eas-c22 VALUE       EAS Y intercept (default: 1.05)
  --h2 VALUE            Legacy common h2 override for all four components
  --h2-outcome VALUE    Legacy common Y1/Y2 h2 override
  --h2-eur-y VALUE      EUR Y background h2 and ALPHA1 scale (default: 0.2)
  --h2-eur-x VALUE      EUR X background h2 and GAMMA1 scale (default: 0.2)
  --h2-eas-x VALUE      EAS X background h2 and GAMMA2 scale (default: 0.2)
  --h2-eas-y VALUE      EAS Y background h2 and ALPHA2 scale (default: 0.2)
  --foreground-gamma-variance-ratio VALUE
                       Active foreground gamma / background per-SNP variance (default: 10)
  --foreground-alpha-variance-ratio VALUE
                       Foreground alpha / background per-SNP variance (default: 1)
  --background-only-fraction VALUE
                       Background-only SNP fraction per chromosome (default: 0.60)
  --background-ld-r2 VALUE
                       Background-only must have raw r2 < VALUE with every foreground
                       SNP in BOTH EUR/EAS within 1 Mb (0 < VALUE <= 1; default: 0.001)
                       Alias: --background-exclusion-r2; env: BACKGROUND_LD_R2
  --paired-effect-streams
                       Pair Gaussian draws by original SNP across background fractions
  --replicates N        Number of replicates (default: 50)
  --thresholds TEXT     Comma-separated P thresholds (default: 5e-5,5e-8)
  --seed N              Master random seed for the entire study (default: 1234)
  --data-root DIR       Existing ViLMA data root
  --work-root DIR       Simulation output root (otherwise nondefault core count / background r2 adds suffixes)
  --python EXE          Python executable
  --no-fema-qc          Skip XMR heterogeneity/BH filtering (default)
  --fema-qc             Enable XMR heterogeneity/BH filtering
  --xmr-priority MODE   FEMA_P (default), EUR_P, EAS_P, both (FEMA/EUR), or all; LD panel is EAS
  --overwrite
  -h, --help

Examples:
  bash run_simulation_study.sh --stage all
  bash run_simulation_study.sh --stage all --method XMR
  bash run_simulation_study.sh --stage methods --method XMR --xmr-priority both
  bash run_simulation_study.sh --stage data --beta 0 0.2 --rho 0.4 --n-aux 40000
  bash run_simulation_study.sh --stage all --n-shared 10
  bash run_simulation_study.sh --stage all --n-shared 50 --n-target-only 5 --n-aux-only 5
  bash run_simulation_study.sh --stage all --background-ld-r2 0.2
  bash run_simulation_study.sh --stage methods --method MRBEE \
    --beta 0,0.2 --rho 0.4,0.8 \
    --n-aux 40000,60000 --n-target 8000,10000
  bash run_simulation_study.sh --stage methods --method XMR --fema-qc
EOF
}

split_grid_values() {
  GRID_VALUES=()
  local token
  local -a pieces
  for token in "$@"; do
    IFS=',' read -r -a pieces <<< "${token}"
    GRID_VALUES+=("${pieces[@]}")
  done
}

export XMR_BG_EURX_EASY_CORR="${XMR_BG_EURX_EASY_CORR:-0.3}"
export XMR_BG_EURY_EASX_CORR="${XMR_BG_EURY_EASX_CORR:-0.1}"
export XMR_EAS_BACKGROUND_CORR="${XMR_EAS_BACKGROUND_CORR:-0.4}"
export XMR_EAS_C11="${XMR_EAS_C11:-1.10}"
export XMR_EAS_C12="${XMR_EAS_C12:-0.15}"
export XMR_EAS_C22="${XMR_EAS_C22:-1.05}"
while (($#)); do
  case "$1" in
    --stage)
      STAGE="$2"; shift 2 ;;
    --method)
      METHOD="$2"; shift 2 ;;
    --beta|--betas)
      BETAS=(); shift
      values=()
      while (($#)) && [[ "$1" != --* ]]; do values+=("$1"); shift; done
      split_grid_values "${values[@]}"; BETAS=("${GRID_VALUES[@]}") ;;
    --background-rho) echo "--background-rho is no longer supported: background rho follows --rhos." >&2; exit 2 ;;
    --background-eurx-easy-correlation) XMR_BG_EURX_EASY_CORR="$2"; shift 2 ;;
    --background-eury-easx-correlation) XMR_BG_EURY_EASX_CORR="$2"; shift 2 ;;
    --eas-background-correlation) XMR_EAS_BACKGROUND_CORR="$2"; shift 2 ;;
    --eas-c11) XMR_EAS_C11="$2"; shift 2 ;;
    --eas-c12) XMR_EAS_C12="$2"; shift 2 ;;
    --eas-c22) XMR_EAS_C22="$2"; shift 2 ;;
    --rho|--rhos)
      RHOS=(); shift
      values=()
      while (($#)) && [[ "$1" != --* ]]; do values+=("$1"); shift; done
      split_grid_values "${values[@]}"; RHOS=("${GRID_VALUES[@]}") ;;
    --n-aux)
      N_AUX_VALUES=(); shift
      values=()
      while (($#)) && [[ "$1" != --* ]]; do values+=("$1"); shift; done
      split_grid_values "${values[@]}"; N_AUX_VALUES=("${GRID_VALUES[@]}") ;;
    --n-target)
      N_TARGET_VALUES=(); shift
      values=()
      while (($#)) && [[ "$1" != --* ]]; do values+=("$1"); shift; done
      split_grid_values "${values[@]}"; N_TARGET_VALUES=("${GRID_VALUES[@]}") ;;
    --n-shared)
      if (($# < 2)) || [[ "$2" == --* ]]; then
        echo "--n-shared requires one positive integer" >&2
        exit 2
      fi
      N_SHARED="$2"; shift 2 ;;
    --n-target-only|--n-aux-only)
      if (($# < 2)) || [[ "$2" == --* ]]; then
        echo "$1 requires one nonnegative integer" >&2
        exit 2
      fi
      if [[ "$1" == --n-target-only ]]; then
        N_TARGET_ONLY="$2"
      else
        N_AUX_ONLY="$2"
      fi
      shift 2 ;;
    --h2)
      EUR_Y_H2="$2"; EUR_X_H2="$2"; EAS_X_H2="$2"; EAS_Y_H2="$2"; shift 2 ;;
    --h2-outcome)
      EUR_Y_H2="$2"; EAS_Y_H2="$2"; shift 2 ;;
    --h2-eur-y)
      EUR_Y_H2="$2"; shift 2 ;;
    --h2-eur-x)
      EUR_X_H2="$2"; shift 2 ;;
    --h2-eas-x)
      EAS_X_H2="$2"; shift 2 ;;
    --h2-eas-y)
      EAS_Y_H2="$2"; shift 2 ;;
    --foreground-gamma-variance-ratio|--foreground-alpha-variance-ratio)
      if (($# < 2)) || [[ "$2" == --* ]]; then
        echo "$1 requires a positive finite number" >&2
        exit 2
      fi
      if [[ "$1" == --foreground-gamma-variance-ratio ]]; then
        FOREGROUND_GAMMA_VARIANCE_RATIO="$2"
      else
        FOREGROUND_ALPHA_VARIANCE_RATIO="$2"
      fi
      shift 2 ;;
    --replicates)
      REPLICATES="$2"; shift 2 ;;
    --background-only-fraction)
      BACKGROUND_ONLY_FRACTION="$2"; shift 2 ;;
    --background-ld-r2|--background-exclusion-r2)
      if (($# < 2)) || [[ "$2" == --* ]]; then
        echo "$1 requires one finite number in (0, 1]" >&2
        exit 2
      fi
      BACKGROUND_LD_R2="$2"; shift 2 ;;
    --paired-effect-streams)
      PAIRED_EFFECT_STREAMS=1; shift ;;
    --thresholds)
      THRESHOLDS="$2"; shift 2 ;;
    --seed)
      MASTER_SEED="$2"; shift 2 ;;
    --data-root)
      DATA_ROOT="$2"; shift 2 ;;
    --work-root)
      WORK_ROOT="$2"; WORK_ROOT_EXPLICIT=1; shift 2 ;;
    --python)
      PYTHON="$2"; shift 2 ;;
    --no-fema-qc)
      NO_FEMA_QC=1; shift ;;
    --fema-qc)
      NO_FEMA_QC=0; shift ;;
    --xmr-priority)
      XMR_PRIORITY="$2"; shift 2 ;;
    --overwrite)
      OVERWRITE=1; shift ;;
    -h|--help)
      usage; exit 0 ;;
    *)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 2 ;;
  esac
done
# Include background settings in the output directory name.
WORK_ROOT="${WORK_ROOT}_bgRhoLinked_uxvy${XMR_BG_EURX_EASY_CORR}_vyux${XMR_BG_EURY_EASX_CORR}_easCorr${XMR_EAS_BACKGROUND_CORR}_C${XMR_EAS_C11}_${XMR_EAS_C12}_${XMR_EAS_C22}"

if [[ ! "${N_SHARED}" =~ ^[1-9][0-9]*$ ]]; then
  echo "--n-shared requires one positive integer: ${N_SHARED}" >&2
  exit 2
fi
if [[ ! "${N_TARGET_ONLY}" =~ ^(0|[1-9][0-9]*)$ ]]; then
  echo "--n-target-only requires one nonnegative integer: ${N_TARGET_ONLY}" >&2
  exit 2
fi
if [[ ! "${N_AUX_ONLY}" =~ ^(0|[1-9][0-9]*)$ ]]; then
  echo "--n-aux-only requires one nonnegative integer: ${N_AUX_ONLY}" >&2
  exit 2
fi

if ! "${PYTHON}" -c '
import math, sys
try:
    value = float(sys.argv[1])
except ValueError:
    raise SystemExit("--background-ld-r2 must be a finite number in (0, 1]")
if not math.isfinite(value) or not 0 < value <= 1:
    raise SystemExit("--background-ld-r2 must be a finite number in (0, 1]")
' "${BACKGROUND_LD_R2}"; then
  exit 2
fi

if ((${#BETAS[@]} == 0 || ${#RHOS[@]} == 0 ||
     ${#N_AUX_VALUES[@]} == 0 || ${#N_TARGET_VALUES[@]} == 0)); then
  echo "--beta, --rho, --n-aux and --n-target cannot be empty" >&2
  exit 2
fi
IFS=',' read -r -a THRESHOLD_VALUES <<< "${THRESHOLDS}"
if ((${#THRESHOLD_VALUES[@]} == 0)); then
  echo "--thresholds cannot be empty" >&2
  exit 2
fi
if [[ ! "${REPLICATES}" =~ ^[1-9][0-9]*$ ]]; then
  echo "--replicates must be a positive integer" >&2
  exit 2
fi
if [[ ! "${MASTER_SEED}" =~ ^[0-9]+$ ]] ||
   ((MASTER_SEED > 2147483646)); then
  echo "--seed must be an integer between 0 and 2147483646" >&2
  exit 2
fi
for n_aux in "${N_AUX_VALUES[@]}"; do
  if [[ ! "${n_aux}" =~ ^[1-9][0-9]*$ ]]; then
    echo "Every --n-aux value must be a positive integer: ${n_aux}" >&2
    exit 2
  fi
done
for n_target in "${N_TARGET_VALUES[@]}"; do
  if [[ ! "${n_target}" =~ ^[1-9][0-9]*$ ]]; then
    echo "Every --n-target value must be a positive integer: ${n_target}" >&2
    exit 2
  fi
done
"${PYTHON}" -c '
import math, sys
split = sys.argv.index("--")
betas = [float(x) for x in sys.argv[1:split]]
rhos = [float(x) for x in sys.argv[split + 1:]]
if not all(math.isfinite(x) for x in betas + rhos):
    raise SystemExit("beta and rho must be finite")
if not all(-1 < x < 1 for x in rhos):
    raise SystemExit("every rho must be in (-1, 1)")
' "${BETAS[@]}" -- "${RHOS[@]}"
"${PYTHON}" -c '
import math, sys
values = [float(x) for x in sys.argv[1:]]
if not all(math.isfinite(x) and 0 < x < 1 for x in values):
    raise SystemExit("all h2 values must be finite and in (0, 1)")
' "${EUR_Y_H2}" "${EUR_X_H2}" "${EAS_X_H2}" "${EAS_Y_H2}"
"${PYTHON}" -c '
import math, sys
names = ("--foreground-gamma-variance-ratio", "--foreground-alpha-variance-ratio")
for name, text in zip(names, sys.argv[1:]):
    try:
        value = float(text)
    except ValueError:
        raise SystemExit(f"{name} must be a positive finite number: {text}")
    if not math.isfinite(value) or value <= 0:
        raise SystemExit(f"{name} must be a positive finite number: {text}")
' "${FOREGROUND_GAMMA_VARIANCE_RATIO}" "${FOREGROUND_ALPHA_VARIANCE_RATIO}"
"${PYTHON}" -c '
import math, sys
value = float(sys.argv[1])
if not math.isfinite(value) or not 0 < value <= 1:
    raise SystemExit("--background-only-fraction must be finite and in (0, 1]")
' "${BACKGROUND_ONLY_FRACTION}"

case "${STAGE}" in
  step3|data|background|methods|evaluate|all|check) ;;
  *) echo "Unknown --stage: ${STAGE}" >&2; exit 2 ;;
esac

METHOD_LOWER="$(printf '%s' "${METHOD}" | tr '[:upper:]' '[:lower:]')"
case "${METHOD_LOWER}" in
  all) METHODS=(XMR CAMeRa MRBEE TEMR MRAPSS); EVALUATION_METHODS="all" ;;
  xmr) METHODS=(XMR); EVALUATION_METHODS="XMR_FEMA" ;;
  camera) METHODS=(CAMeRa); EVALUATION_METHODS="CAMeRa_FEMA" ;;
  mrbee|mrbee_tl)
    METHODS=(MRBEE)
    EVALUATION_METHODS="MRBEE_SOURCE,MRBEE_TARGET,MRBEE_TL" ;;
  temr) METHODS=(TEMR); EVALUATION_METHODS="TEMR" ;;
  mrapss) METHODS=(MRAPSS); EVALUATION_METHODS="MRAPSS_TARGET" ;;
  *) echo "Unknown --method: ${METHOD}" >&2; exit 2 ;;
esac
XMR_QC_ARGS=(--no-fema-qc)
if ((NO_FEMA_QC == 0)); then
  XMR_QC_ARGS=(--fema-qc)
fi
XMR_PRIORITY="$(printf '%s' "${XMR_PRIORITY}" | tr '[:lower:]' '[:upper:]')"
case "${XMR_PRIORITY}" in
  FEMA_P) ;;
  EUR_P)
    if [[ "${EVALUATION_METHODS}" == XMR_FEMA ]]; then
      EVALUATION_METHODS="XMR_EUR_P"
    fi ;;
  EAS_P)
    if [[ "${EVALUATION_METHODS}" == XMR_FEMA ]]; then
      EVALUATION_METHODS="XMR_EAS_P"
    fi ;;
  BOTH)
    if [[ "${EVALUATION_METHODS}" == XMR_FEMA ]]; then
      EVALUATION_METHODS="XMR_FEMA,XMR_EUR_P"
    fi ;;
  ALL)
    if [[ "${EVALUATION_METHODS}" == XMR_FEMA ]]; then
      EVALUATION_METHODS="XMR_FEMA,XMR_EUR_P,XMR_EAS_P"
    fi ;;
  *) echo "--xmr-priority must be FEMA_P, EUR_P, EAS_P, both, or all" >&2; exit 2 ;;
esac

float_label() {
  "${PYTHON}" -c \
    'import sys; print(format(float(sys.argv[1]), ".12g").replace("-", "m").replace(".", "p").replace("+", ""))' \
    "$1"
}

setting_id() {
  printf 'beta%s_rho%s_naux%s_ntarget%s' \
    "$(float_label "$1")" "$(float_label "$2")" "$3" "$4"
}

summary_glob() {
  local sid
  sid="$(setting_id "$1" "$2" "$3" "$4")"
  printf '%s/summary_statistics/%s/xmr_sumstats_*_h2eury%s_eurx%s_easx%s_easy%s_*core%ss%st%sa_*_fixedIV_rep*.tsv.gz' \
    "${WORK_ROOT}" "${sid}" \
    "$(float_label "${EUR_Y_H2}")" \
    "$(float_label "${EUR_X_H2}")" \
    "$(float_label "${EAS_X_H2}")" \
    "$(float_label "${EAS_Y_H2}")" \
    "${N_SHARED}" "${N_TARGET_ONLY}" "${N_AUX_ONLY}"
}

# Different core architectures must not share background/method caches.
if ((WORK_ROOT_EXPLICIT == 0)) &&
   [[ "${N_SHARED}" != "${DEFAULT_N_SHARED}" ||
      "${N_TARGET_ONLY}" != "${DEFAULT_N_TARGET_ONLY}" ||
      "${N_AUX_ONLY}" != "${DEFAULT_N_AUX_ONLY}" ]]; then
  WORK_ROOT="${WORK_ROOT}_core${N_SHARED}s${N_TARGET_ONLY}t${N_AUX_ONLY}a"
fi
if ((WORK_ROOT_EXPLICIT == 0)) &&
   [[ "$(float_label "${BACKGROUND_LD_R2}")" != "$(float_label "${DEFAULT_BACKGROUND_LD_R2}")" ]]; then
  WORK_ROOT="${WORK_ROOT}_bgR2lt$(float_label "${BACKGROUND_LD_R2}")"
fi

# Check before replacing settings.tsv or starting any downstream cached stage.
# Older generated studies store the same rule in their architecture metadata.
"${PYTHON}" -c '
import json, re, sys
from pathlib import Path
root = Path(sys.argv[1])
requested = {"background_ld_r2": float(sys.argv[2]), "background_ld_kb": int(sys.argv[3])}
record = root / "manifests/background_ld_parameters.json"
requested_core = dict(zip(("shared", "target_only", "aux_only"), map(int, sys.argv[4:7])))
core_record = root / "manifests/core_counts.json"

def check_core(path, counts):
    if counts is None or any(counts.get(key) != value for key, value in requested_core.items()):
        raise SystemExit(
            f"{path}: core counts differ from this run "
            f"(stored {counts}; requested {requested_core}). "
            "Use a new --work-root and rerun data/background/methods; --overwrite cannot mix designs."
        )

def check(path, r2, kb):
    if float(r2) != requested["background_ld_r2"] or int(kb) != requested["background_ld_kb"]:
        raise SystemExit(
            f"{path}: background LD parameters differ from this run "
            f"(stored r2={r2}, kb={kb}; requested r2={sys.argv[2]}, kb={sys.argv[3]}). "
            "Use a new --work-root and rerun data/background/methods; --overwrite cannot mix designs."
        )

def check_cross_population_c(path, metadata):
    # Stored four-trait order Y1, X1, X2, Y2; check both cross-population blocks.
    matrix = metadata.get("C")
    order = metadata.get("C_order", ["Y1", "X1", "X2", "Y2"])
    if (order != ["Y1", "X1", "X2", "Y2"]
            or not isinstance(matrix, list) or len(matrix) != 4
            or any(not isinstance(row, list) or len(row) != 4 for row in matrix)):
        raise SystemExit(f"{path}: cannot verify cross-population C; use a new --work-root")
    if any(matrix[i][j] != 0 or matrix[j][i] != 0
           for i in (0, 1) for j in (2, 3)):
        raise SystemExit(
            f"{path}: all EUR/EAS cross-population C entries must be zero. "
            "Use a new --work-root and rerun data/background/methods; --overwrite cannot mix designs."
        )

if record.is_file():
    previous = json.loads(record.read_text(encoding="utf-8"))
    check(record, previous["background_ld_r2"], previous["background_ld_kb"])
if core_record.is_file():
    check_core(core_record, json.loads(core_record.read_text(encoding="utf-8")))
architectures = sorted((root / "true_effects/architectures").glob("*_architecture.json"))
for path in architectures:
    metadata = json.loads(path.read_text(encoding="utf-8"))
    rule = metadata["background_selection"]["background_only_ld_exclusion"]
    check(path, rule["r2_threshold"], rule["window_kb"])
    check_core(path, metadata.get("foreground_selection", {}).get("core_counts"))
design = root / "true_effects/simulation_design.json"
if not architectures and design.is_file():
    scenario = json.loads(design.read_text(encoding="utf-8")).get("scenario", "")
    match = re.search(r"(?:^|_)b([^_]+)MbR2lt([^_]+)(?:_|$)", scenario)
    if match is None:
        raise SystemExit(f"{design}: cannot verify background LD parameters; use a new --work-root")
    def number(label):
        return float(label.replace("m", "-").replace("p", "."))
    check(design, number(match.group(2)), round(1000 * number(match.group(1))))
    core_match = re.search(r"(?:^|_)core([0-9]+)s([0-9]+)t([0-9]+)a(?:_|$)", scenario)
    if core_match is None:
        raise SystemExit(f"{design}: cannot verify core counts; use a new --work-root")
    check_core(design, dict(zip(("shared", "target_only", "aux_only"), map(int, core_match.groups()))))
if design.is_file():
    check_cross_population_c(design, json.loads(design.read_text(encoding="utf-8")))
for path in sorted((root / "summary_statistics").glob("*/xmr_sumstats_*.json")):
    metadata = json.loads(path.read_text(encoding="utf-8"))
    rule = metadata["foreground_design"]
    check(path, rule["background_exclusion_r2"], rule["background_exclusion_kb"])
    check_core(path, {key: rule.get("core_" + key) for key in requested_core})
    check_cross_population_c(path, metadata)
record.parent.mkdir(parents=True, exist_ok=True)
temporary = record.with_suffix(".json.tmp")
temporary.write_text(json.dumps(requested, indent=2) + "\n", encoding="utf-8")
temporary.replace(record)
from tempfile import NamedTemporaryFile
with NamedTemporaryFile(mode="w", encoding="utf-8", dir=core_record.parent,
                        prefix=core_record.name + ".", suffix=".tmp", delete=False) as handle:
    temporary_core = Path(handle.name)
    json.dump(requested_core, handle, indent=2)
    handle.write("\n")
temporary_core.replace(core_record)
' "${WORK_ROOT}" "${BACKGROUND_LD_R2}" "${BACKGROUND_LD_KB}" \
  "${N_SHARED}" "${N_TARGET_ONLY}" "${N_AUX_ONLY}"

echo "Core SNP counts before LD expansion: shared=${N_SHARED}, target-only=${N_TARGET_ONLY}, auxiliary-only=${N_AUX_ONLY}"
echo "Intercept truth: all EUR_Y/EUR_X versus EAS_X/EAS_Y C entries are zero"
echo "Background-only LD rule: raw r2 < ${BACKGROUND_LD_R2} in BOTH EUR/EAS within ${BACKGROUND_LD_KB} kb of every foreground SNP"

mkdir -p \
  "${WORK_ROOT}/manifests" \
  "${WORK_ROOT}/true_effects" \
  "${WORK_ROOT}/summary_statistics" \
  "${WORK_ROOT}/xmr_background" \
  "${WORK_ROOT}/mrbee_aux_background" \
  "${WORK_ROOT}/method_results" \
  "${WORK_ROOT}/evaluation" \
  "${WORK_ROOT}/figures" \
  "${WORK_ROOT}/logs"
export PYTHONPYCACHEPREFIX="${PYTHONPYCACHEPREFIX:-${WORK_ROOT}/.pycache}"

SETTINGS_FILE="${WORK_ROOT}/manifests/settings.tsv"
SETTINGS_TMP="${SETTINGS_FILE}.tmp"
printf 'setting_id\tbeta_true\trho\tn_aux\tn_target\th2_eur_y\th2_eur_x\th2_eas_x\th2_eas_y\tmaster_seed\texpected_replicates\tthresholds\tn_thresholds\tsummary_glob\tbackground_dir\n' > "${SETTINGS_TMP}"
for beta in "${BETAS[@]}"; do
  for rho in "${RHOS[@]}"; do
    for n_aux in "${N_AUX_VALUES[@]}"; do
      for n_target in "${N_TARGET_VALUES[@]}"; do
        sid="$(setting_id "${beta}" "${rho}" "${n_aux}" "${n_target}")"
        glob="$(summary_glob "${beta}" "${rho}" "${n_aux}" "${n_target}")"
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
          "${sid}" "${beta}" "${rho}" "${n_aux}" "${n_target}" \
          "${EUR_Y_H2}" "${EUR_X_H2}" "${EAS_X_H2}" "${EAS_Y_H2}" \
          "${MASTER_SEED}" "${REPLICATES}" "${THRESHOLDS}" "${#THRESHOLD_VALUES[@]}" \
          "${glob}" "${WORK_ROOT}/xmr_background/${sid}" \
          >> "${SETTINGS_TMP}"
      done
    done
  done
done
mv "${SETTINGS_TMP}" "${SETTINGS_FILE}"

# Guard empty-array expansions below for Bash 3.2 with set -u.
OVERWRITE_ARGS=()
if ((OVERWRITE)); then OVERWRITE_ARGS+=(--overwrite); fi
PAIRED_EFFECT_ARGS=()
if ((PAIRED_EFFECT_STREAMS)); then PAIRED_EFFECT_ARGS+=(--paired-effect-streams); fi

run_data() {
  "${PYTHON}" "${SCRIPT_DIR}/step4_generate_xmr_true_effects.py" \
    --ld-dir "${DATA_ROOT}/vilma_ld" \
    --out-dir "${WORK_ROOT}/true_effects" \
    --n-independent-foreground "${N_SHARED}" \
    --n-target-only "${N_TARGET_ONLY}" \
    --n-aux-only "${N_AUX_ONLY}" \
    --n-shared "${N_SHARED}" \
    --background-only-fraction "${BACKGROUND_ONLY_FRACTION}" \
    --h2-eur-y "${EUR_Y_H2}" \
    --h2-eur-x "${EUR_X_H2}" \
    --h2-eas-x "${EAS_X_H2}" \
    --h2-eas-y "${EAS_Y_H2}" \
    --foreground-gamma-variance-ratio "${FOREGROUND_GAMMA_VARIANCE_RATIO}" \
    --foreground-alpha-variance-ratio "${FOREGROUND_ALPHA_VARIANCE_RATIO}" \
    --ldscore-quantiles "${LDSCORE_QUANTILES[@]}" \
    --core-r2 "${CORE_LD_R2}" --core-kb "${CORE_LD_KB}" \
    --ld-r2 "${FOREGROUND_LD_R2}" --ld-kb "${FOREGROUND_LD_KB}" \
    --background-exclusion-r2 "${BACKGROUND_LD_R2}" \
    --background-exclusion-kb "${BACKGROUND_LD_KB}" \
    --betas "${BETAS[@]}" \
    --rhos "${RHOS[@]}" \
    --n-aux "${N_AUX_VALUES[@]}" \
    --n-target "${N_TARGET_VALUES[@]}" \
    --replicates "${REPLICATES}" \
    --seed "${MASTER_SEED}" \
    ${PAIRED_EFFECT_ARGS[@]+"${PAIRED_EFFECT_ARGS[@]}"} \
    ${OVERWRITE_ARGS[@]+"${OVERWRITE_ARGS[@]}"}

  "${PYTHON}" "${SCRIPT_DIR}/step5_generate_xmr_summary_statistics.py" \
    --ld-dir "${DATA_ROOT}/vilma_ld" \
    --truth-glob "${WORK_ROOT}/true_effects/*/xmr_truth_*_h2eury$(float_label "${EUR_Y_H2}")_eurx$(float_label "${EUR_X_H2}")_easx$(float_label "${EAS_X_H2}")_easy$(float_label "${EAS_Y_H2}")_*core${N_SHARED}s${N_TARGET_ONLY}t${N_AUX_ONLY}a_*_fixedIV_rep*.tsv.gz" \
    --out-dir "${WORK_ROOT}/summary_statistics" \
    --n-anchors "${N_SHARED}" \
    --n-target-only "${N_TARGET_ONLY}" \
    --n-aux-only "${N_AUX_ONLY}" \
    --n-shared "${N_SHARED}" \
    --background-only-fraction "${BACKGROUND_ONLY_FRACTION}" \
    --h2-eur-y "${EUR_Y_H2}" \
    --h2-eur-x "${EUR_X_H2}" \
    --h2-eas-x "${EAS_X_H2}" \
    --h2-eas-y "${EAS_Y_H2}" \
    --foreground-gamma-variance-ratio "${FOREGROUND_GAMMA_VARIANCE_RATIO}" \
    --foreground-alpha-variance-ratio "${FOREGROUND_ALPHA_VARIANCE_RATIO}" \
    --ldscore-quantiles "${LDSCORE_QUANTILES[@]}" \
    --core-r2 "${CORE_LD_R2}" --core-kb "${CORE_LD_KB}" \
    --ld-r2 "${FOREGROUND_LD_R2}" --ld-kb "${FOREGROUND_LD_KB}" \
    --background-exclusion-r2 "${BACKGROUND_LD_R2}" \
    --background-exclusion-kb "${BACKGROUND_LD_KB}" \
    --betas "${BETAS[@]}" \
    --rhos "${RHOS[@]}" \
    --n-aux "${N_AUX_VALUES[@]}" \
    --n-target "${N_TARGET_VALUES[@]}" \
    --expected-replicates "${REPLICATES}" \
    --seed "${MASTER_SEED}" \
    ${OVERWRITE_ARGS[@]+"${OVERWRITE_ARGS[@]}"}
  "${PYTHON}" "${SCRIPT_DIR}/build_study_manifests.py" --root "${WORK_ROOT}"
}

run_step3() {
  local args=(
    --ld-dir "${DATA_ROOT}/vilma_ld"
    --chromosomes 2 19
    --progress-every 20
  )
  args+=(${OVERWRITE_ARGS[@]+"${OVERWRITE_ARGS[@]}"})
  "${PYTHON}" "${SCRIPT_DIR}/step3_clump_vilma_snps.py" "${args[@]}"
}

run_background() {
  Rscript --vanilla "${SCRIPT_DIR}/step6_estimate_evaluate_xmr_background.R" \
    --root "${WORK_ROOT}" \
    --betas "${BETAS[@]}" \
    --rhos "${RHOS[@]}" \
    --n-aux "${N_AUX_VALUES[@]}" \
    --n-target "${N_TARGET_VALUES[@]}" \
    --n-anchors "${N_SHARED}" \
    --h2-eur-y "${EUR_Y_H2}" \
    --h2-eur-x "${EUR_X_H2}" \
    --h2-eas-x "${EAS_X_H2}" \
    --h2-eas-y "${EAS_Y_H2}" \
    --expected-replicates "${REPLICATES}" \
    --seed "${MASTER_SEED}" \
    ${OVERWRITE_ARGS[@]+"${OVERWRITE_ARGS[@]}"}
  "${PYTHON}" "${SCRIPT_DIR}/build_study_manifests.py" --root "${WORK_ROOT}"
}

run_one_method() {
  local method="$1" beta="$2" rho="$3" n_aux="$4" n_target="$5"
  local sid glob background script output log common
  sid="$(setting_id "${beta}" "${rho}" "${n_aux}" "${n_target}")"
  glob="$(summary_glob "${beta}" "${rho}" "${n_aux}" "${n_target}")"
  background="${WORK_ROOT}/xmr_background/${sid}"
  output="${WORK_ROOT}/method_results/${sid}/${method}"
  log="${WORK_ROOT}/logs/${sid}_${method}.log"

  shopt -s nullglob
  # Quote the directory so ordinary spaces in --work-root survive globbing.
  local pattern="${glob##*/}"
  files=( "${glob%/*}"/${pattern} )
  shopt -u nullglob
  if ((${#files[@]} != REPLICATES)); then
    echo "${sid}: expected ${REPLICATES} summary files, found ${#files[@]}" >&2
    echo "Glob: ${glob}" >&2
    exit 1
  fi
  common="${files[0]}"
  if [[ ! -f "${background}/rep001/background_parameters.rds" ]]; then
    echo "Missing XMR background parameters: ${background}" >&2
    exit 1
  fi

  case "${method}" in
    XMR) script="run_simulation_XMR.R" ;;
    CAMeRa) script="run_simulation_CAMeRa.R" ;;
    MRBEE) script="run_simulation_MRBEE_TL.R" ;;
    TEMR) script="run_simulation_TEMR.R" ;;
    MRAPSS) script="run_simulation_MRAPSS.R" ;;
  esac

  args=(
    --root "${WORK_ROOT}"
    --data-root "${DATA_ROOT}"
    --setting-id "${sid}"
    --input-glob "${glob}"
    --background-dir "${background}"
    --aux-background-dir "${WORK_ROOT}/mrbee_aux_background/${sid}"
    --output-dir "${output}"
    --common-snps "${common}"
    --beta-true "${beta}"
    --rho "${rho}"
    --n-aux "${n_aux}"
    --n-target "${n_target}"
    --expected-replicates "${REPLICATES}"
    --thresholds "${THRESHOLDS}"
    --python "${PYTHON}"
    --seed "${MASTER_SEED}"
  )
  args+=(${OVERWRITE_ARGS[@]+"${OVERWRITE_ARGS[@]}"})
  if [[ "${method}" == XMR ]]; then
    args+=("${XMR_QC_ARGS[@]}")
    args+=(--xmr-priority "${XMR_PRIORITY}")
  fi

  mkdir -p "${output}"
  echo "Running ${method}: ${sid}"
  Rscript --vanilla "${SCRIPT_DIR}/${script}" "${args[@]}" \
    2>&1 | tee "${log}"
}

run_methods() {
  local beta rho n_aux n_target method
  for beta in "${BETAS[@]}"; do
    for rho in "${RHOS[@]}"; do
      for n_aux in "${N_AUX_VALUES[@]}"; do
        for n_target in "${N_TARGET_VALUES[@]}"; do
          for method in "${METHODS[@]}"; do
            run_one_method \
              "${method}" "${beta}" "${rho}" "${n_aux}" "${n_target}"
          done
        done
      done
    done
  done
}

run_evaluation() {
  Rscript --vanilla "${SCRIPT_DIR}/../../plot/reduced_confounding/evaluate_simulation_study.R" \
    --root "${WORK_ROOT}" \
    --settings "${SETTINGS_FILE}" \
    --methods "${EVALUATION_METHODS}" \
    --xmr-priority "${XMR_PRIORITY}" \
    "${XMR_QC_ARGS[@]}"
}

run_check() {
  for script in \
    step2_compute_vilma_ldscores.py \
    step3_clump_vilma_snps.py \
    build_study_manifests.py \
    step4_generate_xmr_true_effects.py \
    step5_generate_xmr_summary_statistics.py \
    step6_estimate_evaluate_xmr_background.R \
    simulation_method_utils.R \
    run_simulation_XMR.R run_simulation_CAMeRa.R \
    run_simulation_MRBEE_TL.R run_simulation_TEMR.R \
    run_simulation_MRAPSS.R \
    vilma_ld_clump.py; do
    [[ -f "${SCRIPT_DIR}/${script}" ]] || {
      echo "Missing script: ${SCRIPT_DIR}/${script}" >&2; exit 1;
    }
  done
  "${PYTHON}" -m py_compile \
    "${SCRIPT_DIR}/step2_compute_vilma_ldscores.py" \
    "${SCRIPT_DIR}/step3_clump_vilma_snps.py" \
    "${SCRIPT_DIR}/build_study_manifests.py" \
    "${SCRIPT_DIR}/step4_generate_xmr_true_effects.py" \
    "${SCRIPT_DIR}/step5_generate_xmr_summary_statistics.py" \
    "${SCRIPT_DIR}/vilma_ld_clump.py"
  for script in "${SCRIPT_DIR}"/*.R; do
    Rscript --vanilla -e "invisible(parse(file=commandArgs(TRUE)[1]))" "${script}"
  done
  echo "Syntax checks passed. Settings: ${SETTINGS_FILE}"
}

case "${STAGE}" in
  step3) run_step3 ;;
  data) run_data ;;
  background) run_background ;;
  methods) run_methods ;;
  evaluate) run_evaluation ;;
  check) run_check ;;
  all)
    run_data
    run_background
    run_methods
    run_evaluation
    ;;
esac

echo "Completed stage=${STAGE}; work root: ${WORK_ROOT}"
