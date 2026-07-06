#!/usr/bin/env bash
set -euo pipefail

# pymol_load_candidates.sh
#
# Wrapper around pymol_load_candidates.py. Reads a redock_out_sorted.csv
# (from collect_redock_candidates.py) inside DIR and prints/writes a PyMOL
# `load` command per mol_XXXX folder's docked.pdbqt, e.g.:
#
#   load /abs/path/top_candidates/mol_0012/docked.pdbqt, mol_0012
#
# Usage:
#   ./pymol_load_candidates.sh /path/to/top_candidates
#   ./pymol_load_candidates.sh /path/to/top_candidates --out load_candidates.pml
#   ./pymol_load_candidates.sh /path/to/top_candidates | pymol -
#
# All flags are forwarded as-is to pymol_load_candidates.py -- run
# `./pymol_load_candidates.sh --help` for the full flag list (--csv,
# --ligand-file, --folder-col, --out, --quiet).

SCRIPT_DIR="${SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
PY_SCRIPT="${PY_SCRIPT:-${SCRIPT_DIR}/pymol_load_candidates.py}"
PYTHON_BIN="${PYTHON_BIN:-python3}"

if [[ ! -f "${PY_SCRIPT}" ]]; then
  echo "[pymol_load_candidates] error: PY_SCRIPT not found at ${PY_SCRIPT}" >&2
  exit 1
fi

exec "${PYTHON_BIN}" "${PY_SCRIPT}" "$@"
