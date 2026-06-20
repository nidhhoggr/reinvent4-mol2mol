#!/bin/bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

#set -x

#cat   ./results/.steps_done                 # steps completed, vs your TARGET_STEPS
#ls -1 ./results/run_*_1.csv | wc -l         # number of chunks finished
#ls -l ./results/libinvent.chkpt             # checkpoint exists + timestamp
#ls -t ./results/run_*_1.csv | head -3       # most recent chunk files

# is the score actually climbing? mean Score (col 4) per chunk, in order:
#for f in $(ls -v ./results/run_*_1.csv); do
#  awk -F, 'NR>1{s+=$4;n++} END{if(n) printf "%-40s meanScore=%.3f  rows=%d\n", FILENAME, s/n, n}' "$f"
#done

dcol=$(head -1 ${SCRIPT_DIR}/../results/run_00000_1.csv | tr ',' '\n' | grep -n 'docking (raw)' | cut -d: -f1)
for f in $(ls -v ${SCRIPT_DIR}/../results/run_*_1.csv); do
  awk -F, -v d=$dcol 'NR>1{s+=$4; dk+=$d; n++} END{if(n) printf "%-26s meanScore=%.3f  meanDock=%.2f  n=%d\n", FILENAME, s/n, dk/n, n}' "$f"
done
