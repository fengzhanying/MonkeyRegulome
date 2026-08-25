#!/usr/bin/env bash
set -euo pipefail

: "${CCRE_WORKDIR:?Set CCRE_WORKDIR in the config file.}"
: "${SLURM_ARRAY_TASK_ID:?Set SLURM_ARRAY_TASK_ID or run through run_global_ccre.sh.}"

BWAOB="$(command -v bigWigAverageOverBed)"
RDHSBED="${CCRE_WORKDIR}/01.masterlist/monkey_rDHS.bed"
SIGDIR="${CCRE_WORKDIR}/02.signals"
FILELIST="${CCRE_WORKDIR}/logs/bwfile_list.txt"

mkdir -p "${SIGDIR}"

if [[ ! -f "${FILELIST}" ]]; then
  echo "ERROR: ${FILELIST} not found. Run run_global_ccre.sh first." >&2
  exit 1
fi

BW="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${FILELIST}")"
if [[ -z "${BW}" ]]; then
  echo "Task ${SLURM_ARRAY_TASK_ID}: no file at this index, exiting."
  exit 0
fi
if [[ ! -f "${BW}" ]]; then
  echo "ERROR: bigWig file not found: ${BW}" >&2
  exit 1
fi

SAMPLE="$(basename "${BW}" .fc.signal.bw)"
OUTFILE="${SIGDIR}/${SAMPLE}.tab"

if [[ -f "${OUTFILE}" ]]; then
  echo "Already done, skipping: ${SAMPLE}"
  exit 0
fi

echo "Processing [${SLURM_ARRAY_TASK_ID}]: ${SAMPLE}"
"${BWAOB}" "${BW}" "${RDHSBED}" "${OUTFILE}"
echo "Done: ${SAMPLE}"
