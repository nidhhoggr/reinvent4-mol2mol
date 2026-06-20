#!/usr/bin/env bash
# Resilient REINVENT runner.
#
# Runs the RL job in small "chunks" of steps. After each chunk REINVENT writes
# a checkpoint; the next chunk resumes from it. If the process/container is
# killed, just start this script again — it picks up from the last checkpoint.
#
# Use it as the container entrypoint together with a restart policy, e.g.:
#   docker run --gpus all --restart unless-stopped \
#       -v /host/workspace:/workspace reinvent4-gpu \
#       bash /workspace/scripts/run_resilient.sh
#
# Graceful manual stop (lets the current chunk finish its save):
#   docker stop -t 180 <container>
#
# Tunables (override via env): TARGET_STEPS, CHUNK_STEPS
set -uo pipefail

# Activate the env directly (don't rely on ~/.bashrc, which isn't sourced for a
# non-interactive `bash script.sh`, and avoid `conda run` which buffers stdout).
source /opt/conda/etc/profile.d/conda.sh
conda activate reinvent4

RESULTS="/workspace/results"
CFGDIR="/workspace/configs"
PRIOR="$CFGDIR/priors/libinvent.prior"
CHKPT="$RESULTS/libinvent.chkpt"
TEMPLATE="$CFGDIR/libinvent_rl.template.toml"
RUNCFG="$RESULTS/libinvent_rl.run.toml"
PROGRESS="$RESULTS/.steps_done"
DONE="$RESULTS/.complete"

TARGET_STEPS="${TARGET_STEPS:-500}"   # total RL steps you ultimately want
CHUNK_STEPS="${CHUNK_STEPS:-10}"      # steps per checkpointed chunk

mkdir -p "$RESULTS"

if [ -f "$DONE" ]; then
  echo "Run already marked complete ($DONE exists). Idling so the restart"
  echo "policy doesn't loop; run 'docker compose down' to remove the container."
  exec sleep infinity
fi

done_steps=0
[ -f "$PROGRESS" ] && done_steps="$(cat "$PROGRESS")"

echo "Resuming at $done_steps / $TARGET_STEPS steps (chunk size $CHUNK_STEPS)."

while [ "$done_steps" -lt "$TARGET_STEPS" ]; do
  if [ -f "$CHKPT" ]; then
    AGENT="$CHKPT"; USE="true"
  else
    AGENT="$PRIOR"; USE="false"
  fi

  # Per-chunk CSV so REINVENT's truncating open ("w+") never clobbers history.
  CSVPREFIX="$RESULTS/run_$(printf '%05d' "$done_steps")"

  sed -e "s|__AGENT__|$AGENT|" \
      -e "s|__USE__|$USE|" \
      -e "s|__CSVPREFIX__|$CSVPREFIX|" \
      -e "s|__MAXSTEPS__|$CHUNK_STEPS|" \
      "$TEMPLATE" > "$RUNCFG"

  echo ">>> chunk: steps $done_steps -> $((done_steps + CHUNK_STEPS))  (agent=$AGENT)"

  if python -m reinvent "$RUNCFG"; then
    done_steps=$((done_steps + CHUNK_STEPS))
    echo "$done_steps" > "$PROGRESS"
  else
    rc=$?
    echo "Chunk exited non-zero ($rc) — likely interrupted. State is in $CHKPT;"
    echo "restarting this script will resume from there."
    exit "$rc"
  fi
done

touch "$DONE"
echo "Completed $TARGET_STEPS steps. Checkpoint: $CHKPT"
echo "Per-chunk CSVs: $RESULTS/run_*_1.csv  (concatenate for full history)"
echo "Idling so the restart policy won't loop; run 'docker compose down' to stop."
exec sleep infinity
