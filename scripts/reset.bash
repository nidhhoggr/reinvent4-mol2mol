#!/usr/bin/env bash
# Wipe REINVENT run state so the NEXT launch starts fresh from the prior.
# Run this ONLY when you intend a brand-new run — never on normal restarts,
# which are supposed to resume from the checkpoint.
#
#   Host:       bash reset.sh ./workspace/results
#   Container:  bash reset.sh /workspace/results      (default)
set -euo pipefail

RESULTS="${1:-/workspace/results}"

if [ ! -d "$RESULTS" ]; then
  echo "No such directory: $RESULTS"; exit 1
fi

echo "About to delete run state in: $RESULTS"
echo "  checkpoint, .steps_done, .complete, run_*_1.csv, docked artifacts, rendered config"
read -rp "Proceed? This discards all progress. [y/N] " ans
case "$ans" in
  [yY]|[yY][eE][sS]) ;;
  *) echo "Aborted — nothing deleted."; exit 0 ;;
esac

rm -f "$RESULTS"/.steps_done \
      "$RESULTS"/.complete \
      "$RESULTS"/libinvent.chkpt \
      "$RESULTS"/libinvent_rl.run.toml \
      "$RESULTS"/run_*_1.csv \
      "$RESULTS"/[0-9]*docked_* \
      "$RESULTS"/dockstream.log

echo "Clean. The next launch will start fresh (look for 'Resuming at 0 / N steps'"
echo "and 'Agent read from .../libinvent.prior')."
