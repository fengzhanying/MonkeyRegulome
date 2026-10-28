#!/usr/bin/env bash
set -eo pipefail
export LC_ALL=C
if [[ $# -ne 1 || ! -f "$1" ]]; then
  echo "Usage: bash scripts/run_global_ccre.sh config/config.sh" >&2
  exit 1
fi
source "$1"
if [[ -n "${CCRE_CONDA_SH:-}" ]]; then
  source "${CCRE_CONDA_SH}"
  conda activate "${CCRE_CONDA_ENV:-ccre}"
fi
set -u
: "${CCRE_WORKDIR:?Set CCRE_WORKDIR}"
: "${CCRE_MASTERLIST:?Set CCRE_MASTERLIST}"
: "${CCRE_GTF:?Set CCRE_GTF}"
export CCRE_WORKDIR CCRE_MASTERLIST CCRE_GTF
export CCRE_DISTANCE_MODE="${CCRE_DISTANCE_MODE:-legacy}"
export CCRE_PYTHON
CCRE_PYTHON="$(command -v "${CCRE_PYTHON:-python3}")"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXECUTOR="${CCRE_EXECUTOR:-local}"
[[ "$EXECUTOR" == local || "$EXECUTOR" == slurm ]] || { echo "Invalid executor" >&2; exit 1; }
if [[ "$EXECUTOR" == slurm ]]; then
  command -v sbatch >/dev/null
fi
"${CCRE_PYTHON}" "${SCRIPT_DIR}/prepare_run.py"
"${CCRE_PYTHON}" "${SCRIPT_DIR}/00_prepare_masterlist.py"
"${CCRE_PYTHON}" "${SCRIPT_DIR}/01_generate_tss.py"
TOTAL="$(awk 'END {print NR-1}' "${CCRE_WORKDIR}/logs/experiments.tsv")"
if [[ "$EXECUTOR" == slurm ]]; then
  JID="$(sbatch --parsable --export=ALL \
    --array="1-${TOTAL}%${CCRE_ARRAY_CONCURRENCY:-50}" \
    --job-name=ccre_signal -N 1 --mem="${CCRE_SIGNAL_MEM:-5G}" \
    -t "${CCRE_SIGNAL_TIME:-00:30:00}" \
    -o "${CCRE_WORKDIR}/logs/signal_%A_%a.out" \
    -e "${CCRE_WORKDIR}/logs/signal_%A_%a.err" \
    "${SCRIPT_DIR}/02_extract_signals.sh")"
  JID="${JID%%;*}"
  DOWNSTREAM="$(sbatch --parsable --export=ALL --dependency="afterok:${JID}" \
    --job-name=ccre_global -N 1 --mem="${CCRE_DOWNSTREAM_MEM:-30G}" \
    -t "${CCRE_DOWNSTREAM_TIME:-04:00:00}" \
    -o "${CCRE_WORKDIR}/logs/global_%j.out" \
    -e "${CCRE_WORKDIR}/logs/global_%j.err" \
    "${SCRIPT_DIR}/run_downstream.sh")"
  printf "Submitted signal job %s and downstream job %s. Results are pending.\n" "$JID" "$DOWNSTREAM"
else
  for ((IDX=1; IDX<=TOTAL; IDX++)); do
    "${CCRE_PYTHON}" "${SCRIPT_DIR}/extract_signal.py" --index "$IDX"
  done
  bash "${SCRIPT_DIR}/run_downstream.sh"
fi
