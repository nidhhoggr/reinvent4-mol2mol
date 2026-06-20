#!/bin/bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NUM=${1:-100}
echo "Agent,Prior,Target,Score,SMILES,SMILES_state,Input_Scaffold,R-groups,docking,docking (raw),QED,QED (raw),SAScore,SAScore (raw),MW,MW (raw),step" > ${SCRIPT_DIR}/../results/top_candidates.csv
cat ${SCRIPT_DIR}/../results/run_*.csv | sort -t',' -k4 -nr | head -${NUM} >> ${SCRIPT_DIR}/../results/top_candidates.csv
