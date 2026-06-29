#!/bin/sh
# Periodic one-way backup of results/ to Google Drive via rclone.
# Runs INSIDE the rclone/rclone container; loops until the container is stopped,
# then does one final sync on the way out so you don't lose the last interval.
#
# All knobs come from the environment (set them in .env -> compose):
#   RCLONE_REMOTE   rclone remote name from rclone.conf        (default: gdrive)
#   RCLONE_DEST     destination path/folder on the remote      (default: reinvent4-results)
#   SYNC_INTERVAL   seconds between syncs                       (default: 300)
#   BWLIMIT         upload cap, e.g. 2M; "off" = unlimited      (default: off)
#   RCLONE_MODE     copy (never deletes remote) | sync (mirror) (default: copy)
#   RCLONE_EXCLUDES one --exclude pattern                       (default: mol2mol.chkpt)
#   SRC             source dir inside the container             (default: /data/results)
set -u

SRC="${RCLONE_SRC:-results}"
REMOTE="${RCLONE_REMOTE:-gdrive}"
DEST="${RCLONE_DEST:-rclone}"
INTERVAL="${SYNC_INTERVAL:-300}"
BWLIMIT="${BWLIMIT:-off}"
MODE="${RCLONE_MODE:-copy}"
RCLONE_CONFIG="${RCLONE_CONFIG:-configs/rclone.conf}"
# Skip the live rolling checkpoint: it's overwritten in place and a sync can
# catch it mid-write. The immutable checkpoints/step_*.chkpt snapshots (written
# once, never touched again) carry the full history, so excluding the live
# pointer loses nothing and avoids uploading a half-written file.
EXCLUDES="${RCLONE_EXCLUDES:-mol2mol.chkpt}"

run_sync() {
  echo "[rclone] $(date -u +%FT%TZ) $MODE $SRC -> $REMOTE:$DEST (bwlimit=$BWLIMIT)"
  rclone "$MODE" "$SRC" "$REMOTE:$DEST" \
    --exclude "$EXCLUDES" \
    --transfers 2 --checkers 4 --fast-list \
    --bwlimit "$BWLIMIT" \
    --stats 1m --stats-one-line --log-level INFO \
    --config "$RCLONE_CONFIG" \
    || echo "[rclone] returned non-zero (transient?) — will retry next cycle"
}

# Graceful stop: flush once more, then exit cleanly.
trap 'echo "[rclone] stop signal -> final sync"; run_sync; exit 0' TERM INT

echo "[rclone] loop start: every ${INTERVAL}s, mode=$MODE, excluding '$EXCLUDES'"
while true; do
  run_sync
  # backgrounded sleep + wait so the TERM/INT trap fires immediately on stop
  sleep "$INTERVAL" &
  wait $!
done
