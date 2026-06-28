#!/bin/bash
set -eu

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS_DIR="${SCRIPT_DIR}/../results"
OUT="${RESULTS_DIR}/top_candidates.csv"

NUM=100

usage() {
    cat >&2 <<EOF
Usage: $(basename "$0") [-n NUM] <csv-name-or-glob> [more globs...]

  -n NUM   number of top candidates to keep (default: 100)
  -h       show this help

The CSV argument is a filename or glob. Bare names/globs (no slash) are
resolved relative to ${RESULTS_DIR}; arguments containing a slash are used
as-is. All matched files must share an identical header or the script errors.
EOF
    exit "${1:-1}"
}

while getopts ":n:h" opt; do
    case "$opt" in
        n)  NUM="$OPTARG" ;;
        h)  usage 0 ;;
        :)  echo "Error: -$OPTARG requires an argument" >&2; usage ;;
        \?) echo "Error: invalid option -$OPTARG" >&2; usage ;;
    esac
done
shift $((OPTIND - 1))

[[ $# -ge 1 ]] || { echo "Error: no CSV file or glob given" >&2; usage; }
[[ "$NUM" =~ ^[0-9]+$ ]] || { echo "Error: -n must be a positive integer (got '$NUM')" >&2; exit 1; }

# Expand each argument into the concrete file list.
shopt -s nullglob
declare -a files=()
for arg in "$@"; do
    if [[ "$arg" == */* ]]; then
        matches=( $arg )                  # path given: glob/use as-is (CWD-relative or absolute)
    else
        matches=( "${RESULTS_DIR}/"$arg )  # bare name/glob: anchor to the results dir
    fi
    for m in "${matches[@]}"; do
        if [[ "$m" -ef "$OUT" ]]; then continue; fi   # never fold the output file back in
        files+=( "$m" )
    done
done
shopt -u nullglob

[[ ${#files[@]} -ge 1 ]] || { echo "Error: no files matched: $*" >&2; exit 1; }

# Extract the header from the first file and verify every other file matches it.
header=""
ref=""
for f in "${files[@]}"; do
    h="$(head -n 1 "$f")"
    if [[ -z "$ref" ]]; then
        header="$h"; ref="$f"
    elif [[ "$h" != "$header" ]]; then
        echo "Error: header mismatch across glob" >&2
        echo "  $ref" >&2
        echo "    $header" >&2
        echo "  $f" >&2
        echo "    $h" >&2
        exit 1
    fi
done

# Write the shared header, then the data rows (headers stripped) sorted top-N by Score (col 4).
{
    printf '%s\n' "$header"
    for f in "${files[@]}"; do
        tail -n +2 "$f"
    done | sort -t',' -k4 -nr | head -n "$NUM"
} > "$OUT"

echo "Merged ${#files[@]} file(s) -> $OUT (top $NUM by column 4)" >&2
