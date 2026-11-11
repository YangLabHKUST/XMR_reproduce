#!/usr/bin/env bash
# Figure 3 settings: independent effect streams, master seed 1234.
FIXED_ARGS=(
  --data-root "${XMR_SIM_DATA_ROOT:?Set XMR_SIM_DATA_ROOT}"
  --method all --thresholds 5e-5,5e-8 --seed 1234
  --h2-eur-y 0.3 --h2-eur-x 0.3 --h2-eas-x 0.3 --h2-eas-y 0.3
  --n-shared 70 --n-target-only 10 --n-aux-only 10
  --foreground-gamma-variance-ratio 20 --foreground-alpha-variance-ratio 1
  --background-only-fraction 0.60 --background-exclusion-r2 0.1
  --no-fema-qc --xmr-priority FEMA_P
)
