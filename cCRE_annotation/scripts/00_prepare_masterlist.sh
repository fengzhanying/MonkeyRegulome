#!/usr/bin/env bash
set -euo pipefail

: "${CCRE_WORKDIR:?Set CCRE_WORKDIR in the config file.}"
: "${CCRE_MASTERLIST:?Set CCRE_MASTERLIST in the config file.}"

OUTDIR="${CCRE_WORKDIR}/01.masterlist"
OUTBED="${OUTDIR}/monkey_rDHS.bed"
mkdir -p "${OUTDIR}"

echo "Input peak count: $(wc -l < "${CCRE_MASTERLIST}")"

FIRST_CHR="$(awk 'NR==1{print $1}' "${CCRE_MASTERLIST}")"
if [[ "${FIRST_CHR}" == chr* ]]; then
  CHR_PREFIX=""
else
  CHR_PREFIX="chr"
fi

awk -v prefix="${CHR_PREFIX}" '
BEGIN { OFS="\t" }
{
  chrom = prefix $1
  if (chrom ~ /_/ || chrom == "chrM" || chrom == "chrMT") next
  print chrom, $2, $3, chrom"_"$2"_"$3
}' "${CCRE_MASTERLIST}" | sort -k1,1 -k2,2n > "${OUTBED}"

echo "rDHS count after filtering: $(wc -l < "${OUTBED}")"
echo "Output: ${OUTBED}"
