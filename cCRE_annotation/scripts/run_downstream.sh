#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
"${CCRE_PYTHON:-python3}" "${SCRIPT_DIR}/03_calculate_zscores.py"
"${CCRE_PYTHON:-python3}" "${SCRIPT_DIR}/04_classify_ccres.py" --distance-mode "${CCRE_DISTANCE_MODE:-legacy}"
echo "Global cCRE classification completed: ${CCRE_WORKDIR}/04.classification"
