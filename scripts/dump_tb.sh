#!/usr/bin/env bash
#
# dump_tb.sh -- convenience wrapper that runs dump_tb.py inside the reinvent
# container to export TensorBoard scalars to CSV. Runs on the HOST; it shells
# into Docker. No GPU needed -- this just parses the tfevents protobufs.
#
# Usage:
#   ./dump_tb.sh                                    # tb_rl -> results/tb_stats.csv
#   ./dump_tb.sh /workspace/tb_tl                   # override logdir, keep default out
#   ./dump_tb.sh /workspace/tb_rl2 /workspace/out   # override logdir AND outdir
#   LOGDIR=/workspace/tb_tl ./dump_tb.sh            # same via env var
#   IMAGE=reinvent4:latest ./dump_tb.sh            # different image
#
# Unlike redock.sh, dump_tb.py takes POSITIONAL args (logdir, outdir), so args
# you pass REPLACE the defaults rather than appending to them:
#   0 args -> LOGDIR + OUTDIR defaults
#   1 arg  -> that logdir + OUTDIR default
#   2+ args-> passed straight through
#
# Output lands in $OUTDIR inside the container, which maps to <host-cwd>/results
# on the host via the /workspace bind mount -- grab tb_stats.csv from there.
#
# Requirements:
#   - dump_tb.py at ./scripts/dump_tb.py (host) -> /workspace/scripts/dump_tb.py
#   - PYTHON must point at an env that has `tensorboard` (the reinvent4 env does,
#     since it wrote the logs). Confirm once with:
#       docker run --rm htvs-pipeline:latest \
#         /opt/conda/envs/reinvent4/bin/python -c \
#         'import tensorboard; print(tensorboard.__version__)'
#
set -euo pipefail
# ---- overridable defaults (export the var before calling, or edit here) ----
IMAGE="${IMAGE:-htvs-pipeline:latest}"
PYTHON="${PYTHON:-/opt/conda/envs/reinvent4/bin/python}"   # env that has tensorboard
SCRIPT="${SCRIPT:-/workspace/scripts/dump_tb.py}"
LOGDIR="${LOGDIR:-/workspace/tb_rl}"
OUTDIR="${OUTDIR:-/workspace/results}"
WORKDIR_HOST="${WORKDIR_HOST:-$PWD}"        # mounted to /workspace

# Positional args override the LOGDIR/OUTDIR defaults (see header).
case "$#" in
  0) ARGS=("$LOGDIR" "$OUTDIR") ;;
  1) ARGS=("$1" "$OUTDIR") ;;
  *) ARGS=("$@") ;;
esac

echo ">> dump_tb via $IMAGE   logdir=${ARGS[0]}  out=${ARGS[1]}" >&2

exec docker run --rm -v "${WORKDIR_HOST}:/workspace" "$IMAGE" \
  "$PYTHON" "$SCRIPT" "${ARGS[@]}"
