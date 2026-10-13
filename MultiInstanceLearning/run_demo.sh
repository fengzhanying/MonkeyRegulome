#!/usr/bin/env bash
# Train organ-specific MIL models, then predict Brain E130 cCRE labels.
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: bash run_demo.sh [--dry-run]

Run from any directory; all data paths are relative to this script.
  Training inputs: CountMatrix/Monkey*Count.txt
  Training outputs: saved_model/<organ>/final_model.pt, run_config.json,
                    cv_metric_summary.tsv, cv_metrics_by_fold.tsv,
                    cv_training_log.tsv, final_training_log.tsv,
                    final_model_selection_summary.tsv and feature metadata.
  Prediction input: Prediction/Brain/MonkeyBrainE130.prediction_input.tsv
  Prediction output: Prediction/Brain/MonkeyBrainE130.prediction_output.tsv
                     (region_id, predicted_label, max_probability, prob_*).

Environment variables:
  PYTHON_BIN  Python executable with numpy, pandas, torch and scikit-learn
              installed (default: python3).
  DEVICE      auto, cpu or cuda (default: auto). CUDA enables --amp;
              CPU runs without mixed precision.

Examples:
  bash run_demo.sh
  DEVICE=cpu bash run_demo.sh
  PYTHON_BIN=/path/to/python3 DEVICE=cuda bash run_demo.sh
  bash run_demo.sh --dry-run

This runs the supplied matrices with the reference training parameters:
5 folds, up to 50 epochs, early-stopping patience 8. It is a full-data demo,
not a shortened training run. Runtime depends on data size and hardware;
elapsed seconds are printed for every step and for the complete workflow.
Existing model and prediction outputs at the paths above are overwritten.
--dry-run checks dependencies and input files and prints commands without
training or writing outputs. No SLURM installation is required.
EOF
}

DRY_RUN=0
case "${1:-}" in
  --help|-h) usage; exit 0 ;;
  --dry-run) DRY_RUN=1; shift ;;
esac
if [[ $# -ne 0 ]]; then
  usage >&2
  exit 2
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python3}"
DEVICE="${DEVICE:-auto}"

# Validate the environment before starting a potentially long training run.
DEVICE="$("${PYTHON_BIN}" - "${DEVICE}" <<'PY'
import sys
import numpy
import pandas
import sklearn
import torch

device = sys.argv[1]
if device not in {'auto', 'cpu', 'cuda'}:
    sys.exit('DEVICE must be auto, cpu or cuda')
if device == 'auto':
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
if device == 'cuda' and not torch.cuda.is_available():
    sys.exit('CUDA is unavailable in this Python environment; use DEVICE=cpu')
print(device)
PY
)"

shopt -s nullglob
COUNT_MATRICES=("${SCRIPT_DIR}"/CountMatrix/Monkey*Count.txt)
if [[ ${#COUNT_MATRICES[@]} -eq 0 ]]; then
  echo "No training matrices found in ${SCRIPT_DIR}/CountMatrix/" >&2
  exit 1
fi
PREDICTION_INPUT="${SCRIPT_DIR}/Prediction/Brain/MonkeyBrainE130.prediction_input.tsv"
PREDICTION_OUTPUT="${SCRIPT_DIR}/Prediction/Brain/MonkeyBrainE130.prediction_output.tsv"
for INPUT in "${SCRIPT_DIR}/MonkeyMIL_Train.py" "${SCRIPT_DIR}/MonkeyMIL_Pred.py" \
             "${SCRIPT_DIR}/CountMatrix/MonkeyBrainCount.txt" "${PREDICTION_INPUT}" \
             "${COUNT_MATRICES[@]}"; do
  if [[ ! -s "${INPUT}" ]]; then
    echo "Required file is missing or empty: ${INPUT}" >&2
    exit 1
  fi
done

DEVICE_ARGS=(--device "${DEVICE}")
if [[ "${DEVICE}" == cuda ]]; then
  DEVICE_ARGS+=(--amp)
fi

run_step() {
  local started=${SECONDS}
  printf '\nCommand:'
  printf ' %q' "$@"
  printf '\n'
  if [[ "${DRY_RUN}" -eq 0 ]]; then
    "$@"
    printf 'Step completed in %s seconds.\n' "$((SECONDS - started))"
  fi
}

WORKFLOW_STARTED=${SECONDS}
printf 'Python: %s\nDevice: %s\nTraining matrices: %s\n' \
  "${PYTHON_BIN}" "${DEVICE}" "${#COUNT_MATRICES[@]}"

for COUNT_MATRIX in "${COUNT_MATRICES[@]}"; do
  TISSUE="${COUNT_MATRIX##*/}"
  TISSUE="${TISSUE#Monkey}"
  TISSUE="${TISSUE%Count.txt}"
  if [[ -z "${TISSUE}" ]]; then
    echo "Cannot infer organ name from ${COUNT_MATRIX}" >&2
    exit 1
  fi
  printf '\nTraining %s\n' "${TISSUE}"
  run_step "${PYTHON_BIN}" "${SCRIPT_DIR}/MonkeyMIL_Train.py" \
    --count-matrix "${COUNT_MATRIX}" \
    --out-dir "${SCRIPT_DIR}/saved_model/${TISSUE}" \
    --tissue "${TISSUE}" \
    --label-order CA-CTCF CA-H3K4me3 PLS dELS pELS non_cCRE \
    --min-observed-markers 3 \
    --n-splits 5 \
    --learning-rate 3e-3 \
    --lr-scheduler warmup-cosine \
    --warmup-epochs 5 \
    --min-learning-rate 1e-5 \
    --val-fraction 0.1 \
    --seed 20260613 \
    --epochs 50 \
    --patience 8 \
    --batch-size 256 \
    --eval-batch-size 512 \
    --context-dropout-prob 0.0 \
    --single-context-prob 0.0 \
    "${DEVICE_ARGS[@]}"
done

printf '\nPredicting Brain E130\n'
run_step "${PYTHON_BIN}" "${SCRIPT_DIR}/MonkeyMIL_Pred.py" \
  --model "${SCRIPT_DIR}/saved_model/Brain/final_model.pt" \
  --features "${PREDICTION_INPUT}" \
  --out "${PREDICTION_OUTPUT}" \
  --context-column context \
  "${DEVICE_ARGS[@]}"

if [[ "${DRY_RUN}" -eq 1 ]]; then
  printf '\nDry run complete; no models or predictions were written.\n'
else
  printf '\nDemo completed in %s seconds.\nModels: %s/saved_model/\nPredictions: %s\n' \
    "$((SECONDS - WORKFLOW_STARTED))" "${SCRIPT_DIR}" "${PREDICTION_OUTPUT}"
fi
