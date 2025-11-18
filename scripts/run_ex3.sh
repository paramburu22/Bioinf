#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "$PROJECT_ROOT"

EMAIL="${ENTREZ_EMAIL:-}"

if [[ $# -gt 0 ]]; then
  EMAIL="$1"
  shift
fi

if [[ -z "$EMAIL" ]]; then
  cat <<'EOF'
Uso: ./scripts/run_ex3.sh <email> [opciones_extras]

Debe especificarse un email válido (o exportarlo como ENTREZ_EMAIL) para
realizar las descargas por Entrez/NCBI.

Ejemplo:
  ./scripts/run_ex3.sh alumno@universidad.edu --blast-xml data/results/ex02/blast_hbb_remote.xml
EOF
  exit 1
fi

python3 "src/exercises/ex03_msa.py" --email "$EMAIL" "$@"

