#!/usr/bin/env bash
# Interval sidecar: refresh results/tb_stats.csv from the live TB events so the
# rclone backup has an up-to-date view to push to Drive.
#
# Decoupled from REINVENT on purpose: there is no per-step hook, and REINVENT's
# TB writer only flushes ~every 120s, so an interval matches the real data
# cadence without coupling to the training process. Runs in the reinvent image
# (needs tensorboard).
#
# Tunables (env): TB_ROOT, OUTDIR, DUMP_INTERVAL
set -uo pipefail
source /opt/conda/etc/profile.d/conda.sh
conda activate reinvent4

TB_ROOT="${TB_ROOT:-/workspace/results}"     # where REINVENT's TB dir lives (tb_logdir under results)
OUTDIR="${OUTDIR:-/workspace/results}"       # tb_stats.csv lands here (synced by rclone)
INTERVAL="${DUMP_INTERVAL:-120}"             # matches SummaryWriter's flush cadence
SCRIPT="${SCRIPT:-/workspace/scripts/dump_tb.py}"

do_dump() {
  python "$SCRIPT" "$TB_ROOT" "$OUTDIR" \
    || echo "[dump_tb] non-zero (no events flushed yet?) — retry next cycle" >&2
}

trap 'echo "[dump_tb] stop -> final dump"; do_dump; exit 0' TERM INT

echo "[dump_tb] loop: every ${INTERVAL}s, root=$TB_ROOT -> $OUTDIR"
while true; do
  do_dump
  sleep "$INTERVAL" & wait $!     # interruptible sleep so the trap fires promptly
done
