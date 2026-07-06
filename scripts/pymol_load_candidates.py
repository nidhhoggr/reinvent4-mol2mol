#!/usr/bin/env python3
"""
pymol_load_candidates.py

Generates PyMOL `load` commands for each docked.pdbqt referenced in a
redock_out_sorted.csv (the output of collect_redock_candidates.py), one
per mol_XXXX folder, e.g.:

    load /abs/path/to/top_candidates/mol_0012/docked.pdbqt, mol_0012

Commands are emitted in the same order as the CSV (already docking-ranked,
best first). Output goes to stdout by default so you can pipe it straight
into PyMOL, and/or to a .pml file with --out.

Usage:
    ./pymol_load_candidates.py /path/to/top_candidates
    ./pymol_load_candidates.py /path/to/top_candidates --out load_candidates.pml
    ./pymol_load_candidates.py /path/to/top_candidates | pymol -

    # then in PyMOL, or appended to the .pml:
    #   pymol /path/to/top_candidates/load_candidates.pml

Assumes:
  - DIR contains a CSV (default name redock_out_sorted.csv, override with
    --csv) with at least "folder" and "docking" columns, as written by
    collect_redock_candidates.py.
  - Each row's folder (e.g. mol_0012) is a subdirectory of DIR containing
    a --ligand-file (default docked.pdbqt).
"""

import argparse
import csv
import sys
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("dir", type=Path, help="Path to the top_candidates directory (contains mol_XXXX subfolders)")
    p.add_argument("--csv", default=None, type=Path,
                   help="CSV to read (default: <dir>/redock_out_sorted.csv)")
    p.add_argument("--ligand-file", default="docked.pdbqt",
                   help="Filename inside each mol_XXXX folder to load (default: docked.pdbqt)")
    p.add_argument("--folder-col", default="folder", help="Folder column name in the CSV")
    p.add_argument("--out", type=Path, default=None,
                   help="Also write the commands to this .pml file (in addition to stdout)")
    p.add_argument("--quiet", action="store_true",
                   help="Suppress stdout output (useful with --out, or to only get warnings)")
    return p.parse_args()


def main():
    args = parse_args()

    top_dir = args.dir.resolve()
    if not top_dir.is_dir():
        sys.exit(f"error: {top_dir} is not a directory")

    csv_path = (args.csv or (top_dir / "redock_out_sorted.csv")).resolve()
    if not csv_path.exists():
        sys.exit(f"error: CSV not found at {csv_path} (pass --csv to override)")

    with csv_path.open(newline="") as f:
        reader = csv.DictReader(f)
        if args.folder_col not in reader.fieldnames:
            sys.exit(f"error: column '{args.folder_col}' not found in {csv_path}. "
                      f"Available columns: {reader.fieldnames}")
        rows = list(reader)

    lines = []
    for row in rows:
        folder = row[args.folder_col]
        ligand_path = top_dir / folder / args.ligand_file
        if not ligand_path.exists():
            print(f"warning: {ligand_path} not found, skipping {folder}", file=sys.stderr)
            continue
        lines.append(f"load {ligand_path}, {folder}")

    output = "\n".join(lines) + ("\n" if lines else "")

    if not args.quiet:
        sys.stdout.write(output)

    if args.out:
        args.out.write_text(output)
        print(f"wrote {len(lines)} load command(s) to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
