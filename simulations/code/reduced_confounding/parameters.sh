#!/usr/bin/env bash
# Figure 3 settings: independent effect streams, master seed 1234.
export XMR_BG_EURX_EASY_CORR=0.1 XMR_BG_EURY_EASX_CORR=0.1
export XMR_EAS_BACKGROUND_CORR=0.2 XMR_EAS_C11=1.05 XMR_EAS_C12=0.05 XMR_EAS_C22=1.05
RESULT_SUFFIX=_bgRhoLinked_uxvy0.1_vyux0.1_easCorr0.2_C1.05_0.05_1.05
FIXED_ARGS=(
  --data-root "${XMR_SIM_DATA_ROOT:?Set XMR_SIM_DATA_ROOT}"
  --method all --thresholds 5e-5,5e-8 --seed 1234
  --h2-eur-y 0.3 --h2-eur-x 0.3 --h2-eas-x 0.3 --h2-eas-y 0.3
  --n-shared 70 --n-target-only 10 --n-aux-only 10
  --foreground-gamma-variance-ratio 20 --foreground-alpha-variance-ratio 1
  --background-only-fraction 0.60 --background-exclusion-r2 0.1
  --no-fema-qc --xmr-priority FEMA_P
  --eas-background-correlation 0.2 --eas-c11 1.05 --eas-c12 0.05 --eas-c22 1.05
  --background-eurx-easy-correlation 0.1 --background-eury-easx-correlation 0.1
)
