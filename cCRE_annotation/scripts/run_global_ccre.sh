#!/usr/bin/env bash
set -euo pipefail

if [[ $# -gt 1 ]]; then
  echo "Usage: bash scripts/run_global_ccre.sh [config/config.sh]" >&2
  exit 1
fi

if [[ $# -eq 1 ]]; then
  # shellcheck disable=SC1090
  source "$1"
fi

: "${CCRE_WORKDIR:?Set CCRE_WORKDIR in the config file.}"
: "${CCRE_MASTERLIST:?Set CCRE_MASTERLIST in the config file.}"
: "${CCRE_GTF:?Set CCRE_GTF in the config file.}"
: "${CCRE_ATAC_BW_DIR:?Set CCRE_ATAC_BW_DIR in the config file.}"
: "${CCRE_CUTTAG_BW_DIR:?Set CCRE_CUTTAG_BW_DIR in the config file.}"

SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FILELIST="${CCRE_WORKDIR}/logs/bwfile_list.txt"

mkdir -p \
  "${CCRE_WORKDIR}/01.masterlist" \
  "${CCRE_WORKDIR}/02.signals" \
  "${CCRE_WORKDIR}/03.zscores" \
  "${CCRE_WORKDIR}/04.classification" \
  "${CCRE_WORKDIR}/logs"

if [[ -n "${CCRE_CONDA_SH:-}" ]]; then
  # shellcheck disable=SC1090
  source "${CCRE_CONDA_SH}"
  conda activate "${CCRE_CONDA_ENV:-ucsctools}"
fi

command -v bigWigAverageOverBed >/dev/null
command -v bedtools >/dev/null
command -v python3 >/dev/null

echo "Step 1: prepare macaque cCRE/rDHS anchor list"
bash "${SCRIPTS_DIR}/00_prepare_masterlist.sh"

echo "Step 2: generate TSS annotations"
python3 "${SCRIPTS_DIR}/01_generate_tss.py"

echo "Step 3: collect ATAC/CUT&Tag bigWig files"
: > "${FILELIST}"
find "${CCRE_ATAC_BW_DIR}" -maxdepth 1 -type f -name 'ATAC-*.fc.signal.bw' | sort >> "${FILELIST}"
find "${CCRE_CUTTAG_BW_DIR}" -maxdepth 1 -type f -name 'CUT-Tag-*-H3K4me3.fc.signal.bw' | sort >> "${FILELIST}"
find "${CCRE_CUTTAG_BW_DIR}" -maxdepth 1 -type f -name 'CUT-Tag-*-H3K27ac.fc.signal.bw' | sort >> "${FILELIST}"
find "${CCRE_CUTTAG_BW_DIR}" -maxdepth 1 -type f -name 'CUT-Tag-*-CTCF.fc.signal.bw' | sort >> "${FILELIST}"

TOTAL="$(wc -l < "${FILELIST}")"
echo "Found ${TOTAL} bigWig files"
if [[ "${TOTAL}" -eq 0 ]]; then
  echo "ERROR: no bigWig files found" >&2
  exit 1
fi

if command -v sbatch >/dev/null; then
  echo "Step 4: submit signal extraction array"
  JID="$(
    sbatch --parsable \
      --array="1-${TOTAL}%${CCRE_ARRAY_CONCURRENCY:-50}" \
      --job-name=ccre_signal \
      -N 1 --mem="${CCRE_SIGNAL_MEM:-5G}" -t "${CCRE_SIGNAL_TIME:-00:30:00}" \
      -o "${CCRE_WORKDIR}/logs/signal_%a.out" \
      -e "${CCRE_WORKDIR}/logs/signal_%a.err" \
      "${SCRIPTS_DIR}/02_extract_signals.sh"
  )"

  echo "Step 5: submit downstream z-score and classification job"
  sbatch \
    --dependency="afterok:${JID}" \
    --job-name=ccre_global \
    -N 1 --mem="${CCRE_DOWNSTREAM_MEM:-30G}" -t "${CCRE_DOWNSTREAM_TIME:-04:00:00}" \
    -o "${CCRE_WORKDIR}/logs/global.out" \
    -e "${CCRE_WORKDIR}/logs/global.err" \
    --wrap="python3 '${SCRIPTS_DIR}/03_calculate_zscores.py' && python3 '${SCRIPTS_DIR}/04_classify_ccres.py'"
else
  echo "Step 4: run signal extraction serially because sbatch is unavailable"
  for IDX in $(seq 1 "${TOTAL}"); do
    SLURM_ARRAY_TASK_ID="${IDX}" bash "${SCRIPTS_DIR}/02_extract_signals.sh"
  done

  echo "Step 5: compute z-scores and classify global cCREs"
  python3 "${SCRIPTS_DIR}/03_calculate_zscores.py"
  python3 "${SCRIPTS_DIR}/04_classify_ccres.py"
fi

echo "Done. Global outputs are in ${CCRE_WORKDIR}/04.classification"
