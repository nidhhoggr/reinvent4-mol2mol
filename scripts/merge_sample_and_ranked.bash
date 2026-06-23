#!/usr/bin/env bash
#
# merge_scores.sh — join ranked docking results with the sampled SMILES table.
#
# redock.csv is sorted by docking score; its `directory` column (mol_NNNN) maps
# positionally into sampled.csv, where:
#
#     line in sampled.csv (header = line 1) = mol_index + 2
#     i.e. mol_0000 -> line 2, mol_0001 -> line 3, ...
#
# Output keeps redock.csv's row order (i.e. ranked by score) and appends the
# sampled.csv columns. Rows with no matching sampled line get blank fields.
#
# Usage:   ./merge_scores.sh <redock.csv> <sampled.csv> [output.csv]
# Default output: merged.csv
#
set -euo pipefail

if [[ $# -lt 2 ]]; then
    echo "Usage: $0 <redock.csv> <sampled.csv> [output.csv]" >&2
    exit 1
fi

REDOCK="$1"
SAMPLED="$2"
OUT="${3:-merged.csv}"

for f in "$REDOCK" "$SAMPLED"; do
    [[ -f "$f" ]] || { echo "ERROR: '$f' not found" >&2; exit 1; }
done

# NOTE: sampled.csv is passed FIRST to awk so it loads before redock is streamed.
awk -F, -v OFS=, '
    # ---- first file: sampled.csv -> index rows by molecule id ----
    NR==FNR {
        if (FNR==1) { shead=$0; scols=NF; next }   # header
        sampled[FNR-2]=$0                          # FNR=2 -> mol_0000
        next
    }
    # ---- second file: redock.csv ----
    FNR==1 { print $0, shead; next }               # combined header
    {
        n=$2; sub(/^mol_/,"",n); n+=0              # mol_0148 -> 148 (numeric)
        if (n in sampled) {
            print $0, sampled[n]
        } else {
            pad=""; for (i=0;i<scols;i++) pad=pad OFS
            print $0 pad                            # blank sampled cols
            nmiss++
        }
    }
    END { if (nmiss>0)
            printf("WARN: %d redock row(s) had no sampled match\n", nmiss) > "/dev/stderr" }
' "$SAMPLED" "$REDOCK" > "$OUT"

echo "Wrote $OUT ($(($(wc -l < "$OUT") - 1)) data rows)." >&2
