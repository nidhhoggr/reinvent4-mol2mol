#!/usr/bin/env bash
#
# rank_scores.sh — collect AutoDock Vina docking scores from mol_* directories
# and write a CSV ranked by best (most negative) score first.
#
# Usage:   ./rank_scores.sh <input_dir> [output.csv]
#   <input_dir>   directory containing the mol_* subdirectories (required)
#   [output.csv]  output file (default: docking_scores.csv in the CWD)
#
# Examples:
#   ./rank_scores.sh redock
#   ./rank_scores.sh /workspace/results/redock my_results.csv
#
set -euo pipefail
shopt -s nullglob

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <input_dir> [output.csv]" >&2
    exit 1
fi

INDIR="${1%/}"                       # strip any trailing slash
OUT="${2:-docking_scores.csv}"
PDBQT="docked.pdbqt"                 # filename inside each mol_* dir to read

if [[ ! -d "$INDIR" ]]; then
    echo "ERROR: '$INDIR' is not a directory" >&2
    exit 1
fi

tmp="$(mktemp)"
trap 'rm -f "$tmp"' EXIT

dirs_seen=0
dirs_ok=0

for dir in "$INDIR"/mol_*/; do
    dir="${dir%/}"                   # strip trailing slash
    name="$(basename "$dir")"        # e.g. mol_2498, not the full path
    dirs_seen=$((dirs_seen + 1))
    file="$dir/$PDBQT"

    if [[ ! -s "$file" ]]; then
        echo "WARN: $file missing or empty — skipped" >&2
        continue
    fi

    # All pose scores in file order. Vina writes the best pose as MODEL 1,
    # so the first value is the best score for that ligand.
    mapfile -t scores < <(grep "VINA RESULT" "$file" | awk '{print $4}')
    if [[ ${#scores[@]} -eq 0 ]]; then
        echo "WARN: no VINA RESULT line in $file — skipped" >&2
        continue
    fi

    best="${scores[0]}"
    second="${scores[1]:-NA}"
    num_poses="${#scores[@]}"

    torsions="$(grep -m1 'active torsions' "$file" | awk '{print $2}')"
    torsions="${torsions:-NA}"

    # Gap to 2nd pose: positive = clear winner over the runner-up
    if [[ "$second" != "NA" ]]; then
        gap="$(awk "BEGIN{printf \"%.3f\", $second - $best}")"
    else
        gap="NA"
    fi

    printf '%s,%s,%s,%s,%s,%s\n' \
        "$name" "$best" "$num_poses" "$torsions" "$second" "$gap" >> "$tmp"
    dirs_ok=$((dirs_ok + 1))
done

if [[ $dirs_seen -eq 0 ]]; then
    echo "ERROR: no mol_* directories found in '$INDIR'" >&2
    exit 1
fi

# Header
echo "rank,directory,best_score,num_poses,num_torsions,second_best_score,gap_to_second" > "$OUT"

# Sort by best_score (field 2). -g = general numeric, so -9.1 sorts before -5.2,
# i.e. the strongest binder (most negative) ends up at the top = rank 1.
# Want literal numeric DESCENDING (-5.2 first) instead? add -r to the sort.
sort -t, -k2,2g "$tmp" | awk -F, 'BEGIN{OFS=","}{print NR, $0}' >> "$OUT"

echo "Wrote $OUT — $dirs_ok of $dirs_seen directories had scores." >&2
