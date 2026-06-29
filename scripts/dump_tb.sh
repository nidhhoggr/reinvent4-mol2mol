#!/usr/bin/env bash
#
# dump_tb.sh -- convenience wrapper that runs dump_tb.py inside the reinvent
# container to export TensorBoard scalars to CSV. Runs on the HOST; it shells
# into Docker. No GPU needed -- this just parses the tfevents protobufs.
#
# dump_tb.py stitches the resilient loop's per-chunk dirs (results/run_<N>_tb)
# into ONE continuous global-step series, writing results/tb_stats.csv. With
# EMIT_TB=1 it also drops results/tb_merged/ -- a single events file with the
# stitched/renumbered steps so you can view one unbroken curve in TensorBoard
# instead of N overlaid per-chunk runs.
#
# Usage:
#   ./dump_tb.sh                                    # results -> results/tb_stats.csv
#   ./dump_tb.sh /workspace/tb_tl                   # override logdir, keep default out
#   ./dump_tb.sh /workspace/results /workspace/out  # override logdir AND outdir
#   LOGDIR=/workspace/tb_tl ./dump_tb.sh            # same via env var
#   EMIT_TB=1 ./dump_tb.sh                          # also write results/tb_merged/ for TensorBoard
#   IMAGE=reinvent4:latest ./dump_tb.sh            # different image
#
# Unlike redock.sh, dump_tb.py takes POSITIONAL args (logdir, outdir), so args
# you pass REPLACE the defaults rather than appending to them:
#   0 args -> LOGDIR + OUTDIR defaults
#   1 arg  -> that logdir + OUTDIR default
#   2+ args-> passed straight through
# The --emit-tb flag is appended separately via EMIT_TB=1 (don't pass it as a
# positional, or it'd be read as the logdir).
#
# Output lands in $OUTDIR inside the container, which maps to <host-cwd>/results
# on the host via the /workspace bind mount -- grab tb_stats.csv (and tb_merged/)
# from there. View the merged curve with:
#   tensorboard --logdir <outdir>/tb_merged --host 0.0.0.0 --port 6006
#
# Requirements:
#   - dump_tb.py at ./scripts/dump_tb.py (host) -> /workspace/scripts/dump_tb.py
#   - PYTHON must point at an env with `tensorboard` (the reinvent4 env, since it
#     wrote the logs). EMIT_TB also needs `torch` (same env). Confirm once with:
#       docker run --rm htvs-pipeline:latest \
#         /opt/conda/envs/reinvent4/bin/python -c \
#         'import tensorboard, torch; print(tensorboard.__version__, torch.__version__)'
#
set -euo pipefail
# ---- overridable defaults (export the var before calling, or edit here) ----
IMAGE="${IMAGE:-htvs-pipeline:latest}"
PYTHON="${PYTHON:-/opt/conda/envs/reinvent4/bin/python}"   # env that has tensorboard (+ torch for EMIT_TB)
SCRIPT="${SCRIPT:-/workspace/scripts/dump_tb.py}"
LOGDIR="${LOGDIR:-/workspace/results}"      # chunk TB dirs (run_*_tb) live under results now
OUTDIR="${OUTDIR:-/workspace/results}"
WORKDIR_HOST="${WORKDIR_HOST:-$PWD}"        # mounted to /workspace

# Positional args override the LOGDIR/OUTDIR defaults (see header).
case "$#" in
  0) ARGS=("$LOGDIR" "$OUTDIR") ;;
  1) ARGS=("$1" "$OUTDIR") ;;
  *) ARGS=("$@") ;;
esac

# EMIT_TB=1 appends the --emit-tb flag (writes <outdir>/tb_merged/).
EXTRA=()
[ "${EMIT_TB:-0}" = "1" ] && EXTRA+=(--emit-tb)

echo ">> dump_tb via $IMAGE   logdir=${ARGS[0]}  out=${ARGS[1]}  emit_tb=${EMIT_TB:-0}" >&2

exec docker run --rm -v "${WORKDIR_HOST}:/workspace" "$IMAGE" \
  "$PYTHON" "$SCRIPT" "${ARGS[@]}" ${EXTRA[@]+"${EXTRA[@]}"}
