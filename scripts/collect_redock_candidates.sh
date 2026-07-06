#!/usr/bin/env bash
set -euo pipefail

# collect_redock_candidates.sh
#
# Wrapper around collect_redock_candidates.py, which post-processes the
# output of vina_redock.py / vina_redock.bash (nidhhoggr/reinvent4-mol2mol):
#
#   https://github.com/nidhhoggr/reinvent4-mol2mol/blob/master/scripts/vina_redock.bash
#
# vina_redock.bash docks a --csv of SMILES and writes, under its OUT dir,
# one mol_{i:04d} folder per input row (i = 0-based row index), each with a
# docked.pdbqt containing Vina's "REMARK VINA RESULT" lines. This script
# expects exactly that layout: point RESULTS_DIR at vina_redock.bash's OUT,
# and CSV at the same file you passed to it (its --csv / positional CSV
# path), so row order lines up with folder index order.
#
# It matches each row's SMILES back to its mol_XXXX folder, pulls the best
# (first) Vina score out of docked.pdbqt, moves folders passing THRESHOLD
# into OUT_DIR, and writes a docking-ranked CSV.
#
# All paths (CSV, RESULTS_DIR, OUT_DIR, OUT_CSV) are relative to WORKDIR
# (defaults to $PWD). collect_redock_candidates.py itself is stdlib-only,
# so no venv/requirements are needed beyond python3.
#
# Usage:
#   # 1) plain CLI flags, forwarded as-is to collect_redock_candidates.py:
#   ./collect_redock_candidates.sh --csv path/to/dock_remaining.csv \
#       --results-dir path/to/pt2 --threshold -11.5 \
#       --out-dir path/to/pt2/top_candidates --out-csv path/to/pt2/redock_out_sorted.csv
#
#   # 2) env-var form (only used when no CLI flags are given):
#   ./collect_redock_candidates.sh
#   THRESHOLD=-11.5 CSV=dock_remaining.csv ./collect_redock_candidates.sh
#   RESULTS_DIR=/workspace/results/redock ./collect_redock_candidates.sh
#   EXTRA_ARGS="--dry-run" ./collect_redock_candidates.sh
#   EXTRA_ARGS="--copy" ./collect_redock_candidates.sh   # copy instead of move

# --- overridable defaults ---------------------------------------------------
SCRIPT_DIR="${SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
PY_SCRIPT="${PY_SCRIPT:-${SCRIPT_DIR}/collect_redock_candidates.py}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
WORKDIR="${WORKDIR:-$(pwd)}"

CSV="${CSV:-dock_remaining.csv}"           # same CSV given to vina_redock.bash
RESULTS_DIR="${RESULTS_DIR:-.}"            # vina_redock.bash's OUT dir
THRESHOLD="${THRESHOLD:--11.0}"
OUT_DIR="${OUT_DIR:-top_candidates}"
OUT_CSV="${OUT_CSV:-redock_out_sorted.csv}"
SMILES_COL="${SMILES_COL:-SMILES}"
EXTRA_ARGS="${EXTRA_ARGS:-}"                # e.g. "--dry-run" or "--copy"
# -----------------------------------------------------------------------------

if [[ ! -f "${PY_SCRIPT}" ]]; then
  echo "[collect_redock_candidates] error: PY_SCRIPT not found at ${PY_SCRIPT}" >&2
  exit 1
fi

cd "${WORKDIR}"

if [[ $# -gt 0 ]]; then
  # CLI flags were passed directly (e.g. --csv ... --results-dir ...):
  # forward them as-is instead of the env-var-built defaults above, which
  # would otherwise silently override/ignore them.
  exec "${PYTHON_BIN}" "${PY_SCRIPT}" "$@"
fi

# shellcheck disable=SC2086
exec "${PYTHON_BIN}" "${PY_SCRIPT}" \
  --csv "${CSV}" \
  --results-dir "${RESULTS_DIR}" \
  --threshold "${THRESHOLD}" \
  --out-dir "${OUT_DIR}" \
  --out-csv "${OUT_CSV}" \
  --smiles-col "${SMILES_COL}" \
  ${EXTRA_ARGS}
