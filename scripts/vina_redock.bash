#!/usr/bin/env bash
#
# redock.sh -- convenience wrapper that runs vina_redock.py inside the
# htvs-pipeline container with all the boilerplate (image, mount, binaries,
# receptor, search box) filled in. Runs on the HOST; it shells into Docker.
#
# Usage:
#   ./redock.sh "Cc1ccc(-c2cccc(CNC3...)n2)cc1"          # one SMILES
#   ./redock.sh --csv /workspace/results/combined.csv     # whole CSV column
#   EXH=16 ./redock.sh "CCO..."                           # override exhaustiveness
#   CENTER="2.0 3.0 18.0" ./redock.sh "CCO..."            # override the box
#   WORKERS=10 ./redock.sh --csv /workspace/results/combined.csv   # 10 ligands at once
#
# Any extra args are passed straight through to vina_redock.py, and because
# they come last, they override the defaults below (argparse takes the last value).
#
# Parallelism: each docking stays single-threaded (--cpu 1, --seed 42) so scores
# stay reproducible against your DockStream run. WORKERS controls how many ligands
# dock concurrently in the --csv path -- set it to your core count (e.g. 10).
# It has no effect when docking a single SMILES.
#
# Requirements:
#   - vina_redock.py placed at ./scripts/vina_redock.py (host) so it appears at
#     /workspace/scripts/vina_redock.py inside the container.
#   - confirm the vina/obabel paths once with:
#       docker run --rm htvs-pipeline:latest bash -lc \
#         'for e in reinvent4 DockStream; do echo "== $e =="; \
#          ls /opt/conda/envs/$e/bin | grep -xE "vina|obabel"; done'
#
set -euo pipefail
# ---- overridable defaults (export the var before calling, or edit here) ----
IMAGE="${IMAGE:-htvs-pipeline:latest}"
RECEPTOR="${RECEPTOR:-/workspace/docking_setup/receptor.pdbqt}"
CENTER="${CENTER:-1.56 2.56 18.4}"          # must match your dockstream_config.json
SIZE="${SIZE:-22 22 22}"
EXH="${EXH:-8}"                             # must match your config's --exhaustiveness
WORKERS="${WORKERS:-1}"                     # parallel ligands for --csv (each still --cpu 1)
VINA="${VINA:-/opt/conda/envs/reinvent4/bin/vina}"
OBABEL="${OBABEL:-/opt/conda/envs/DockStream/bin/obabel}"
PYTHON="${PYTHON:-/opt/conda/envs/DockStream/bin/python}"   # env that has RDKit + obabel
SCRIPT="${SCRIPT:-/workspace/scripts/vina_redock.py}"
OUT="${OUT:-/workspace/results/redock}"
WORKDIR_HOST="${WORKDIR_HOST:-$PWD}"        # mounted to /workspace
if [ "$#" -eq 0 ]; then
  echo "usage: $0 \"<SMILES>\"   |   $0 --csv <path> [--smiles-col COL]" >&2
  exit 1
fi
echo ">> redock via $IMAGE   exh=$EXH  workers=$WORKERS  center=($CENTER)  size=($SIZE)" >&2
# CENTER/SIZE are intentionally unquoted so each splits into three arguments.
exec docker run --rm -v "${WORKDIR_HOST}:/workspace" "$IMAGE" \
  "$PYTHON" "$SCRIPT" \
  --receptor "$RECEPTOR" \
  --center $CENTER \
  --size $SIZE \
  --exhaustiveness "$EXH" \
  --workers "$WORKERS" \
  --vina "$VINA" \
  --obabel "$OBABEL" \
  --out "$OUT" \
  "$@"
