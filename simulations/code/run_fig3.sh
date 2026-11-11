#!/usr/bin/env bash
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=$(cd "$HERE/../.." && pwd)
usage() {
  cat <<'EOF'
Usage:
  bash simulations/code/run_fig3.sh --all [--dry-run]
  bash simulations/code/run_fig3.sh PROFILE TARGET_N [AUX_N] [RHO] [--dry-run]

PROFILE: baseline or reduced_confounding
TARGET_N: 6000 or 8000; AUX_N: 20000, 40000 (default), or 60000
RHO: 0.2, 0.5, or 0.8 (default)
--all runs the configurations in simulations/fig3_settings.tsv.
--dry-run prints commands without generating data or fitting models.
Prepare the LD inputs first, as described in simulations/README.md.
EOF
}
DRY_RUN=0
POSITIONAL=()
for ARG in "$@"; do
  case "$ARG" in
    --dry-run) DRY_RUN=1 ;;
    --help|-h) usage; exit 0 ;;
    *) POSITIONAL+=("$ARG") ;;
  esac
done
set -- "${POSITIONAL[@]}"
if [[ "${1:-}" == --all ]]; then
  [[ $# == 1 ]] || { usage >&2; exit 2; }
  while IFS=$'\t' read -r PROFILE TARGET_N AUX_N RHO; do
    [[ "$PROFILE" == profile ]] && continue
    if [[ "$DRY_RUN" == 1 ]]; then
      bash "$HERE/run_fig3.sh" "$PROFILE" "$TARGET_N" "$AUX_N" "$RHO" --dry-run
    else
      bash "$HERE/run_fig3.sh" "$PROFILE" "$TARGET_N" "$AUX_N" "$RHO"
    fi
  done < "$REPO/simulations/fig3_settings.tsv"
  exit 0
fi
[[ $# -ge 2 && $# -le 4 ]] || { usage >&2; exit 2; }
PROFILE=$1
TARGET_N=$2
AUX_N=${3:-40000}
RHO=${4:-0.8}
[[ "$PROFILE" == baseline || "$PROFILE" == reduced_confounding ]] || { usage >&2; exit 2; }
[[ "$TARGET_N" == 6000 || "$TARGET_N" == 8000 ]] || { usage >&2; exit 2; }
[[ "$AUX_N" == 20000 || "$AUX_N" == 40000 || "$AUX_N" == 60000 ]] || { usage >&2; exit 2; }
[[ "$RHO" == 0.2 || "$RHO" == 0.5 || "$RHO" == 0.8 ]] || { usage >&2; exit 2; }
export XMR_SIM_DATA_ROOT="${XMR_SIM_DATA_ROOT:-$REPO/simulations/data}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1
export R_DATATABLE_NUM_THREADS=1 OMP_THREAD_LIMIT=1 MPLBACKEND=Agg
export GOTO_NUM_THREADS=1 RCPP_PARALLEL_NUM_THREADS=1 TBB_NUM_THREADS=1
unset N_SHARED N_TARGET_ONLY N_AUX_ONLY BACKGROUND_LD_R2
unset XMR_BG_EURX_EASY_CORR XMR_BG_EURY_EASX_CORR
unset XMR_EAS_BACKGROUND_CORR XMR_EAS_C11 XMR_EAS_C12 XMR_EAS_C22
source "$HERE/$PROFILE/parameters.sh"
RHO_LABEL=${RHO/./p}
WORK_ROOT="$REPO/simulations/results/$PROFILE/naux${AUX_N}_ntarget${TARGET_N}_rho${RHO_LABEL}"
for STAGE in data background methods; do
  COMMAND=(bash "$HERE/$PROFILE/run_simulation_study.sh" "${FIXED_ARGS[@]}"
    --stage "$STAGE" --work-root "$WORK_ROOT"
    --beta 0 0.1 0.2 0.3 0.4 --rho "$RHO" --n-aux "$AUX_N"
    --n-target "$TARGET_N" --replicates 50)
  if [[ "$DRY_RUN" == 1 ]]; then
    printf '%q ' "${COMMAND[@]}"; printf '\n'
  else
    "${COMMAND[@]}"
  fi
done
