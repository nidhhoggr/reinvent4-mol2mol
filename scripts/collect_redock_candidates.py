#!/usr/bin/env python3
"""
collect_redock_candidates.py

Post-processes the output of vina_redock.py / vina_redock.bash
(nidhhoggr/reinvent4-mol2mol, scripts/vina_redock.bash):

  https://github.com/nidhhoggr/reinvent4-mol2mol/blob/master/scripts/vina_redock.bash

That script docks a --csv of SMILES and writes, under its --out directory,
one folder per input row named mol_{i:04d} (i = the row's 0-based enumerate
index in the CSV), each containing docked.pdbqt with Vina's
"REMARK VINA RESULT:" lines -- the same file/line format vina_redock.py
itself parses (line.startswith("REMARK VINA RESULT"), best pose first).

This script matches each CSV row's SMILES back to its mol_XXXX folder (by
row order == folder index order, exactly how vina_redock.py created them),
pulls the best score out of docked.pdbqt, moves folders that beat a
threshold into a top_candidates/ dir, and writes a ranked CSV.

Usage:
    ./collect_redock_candidates.py \
        --csv dock_remaining.csv \
        --results-dir . \
        --threshold -11.0 \
        --out-dir top_candidates \
        --out-csv redock_out_sorted.csv

Assumes:
  - dock_remaining.csv has a header row with a "SMILES" column (matching
    --smiles-col, and the same file/column passed to vina_redock.py's --csv).
  - Result folders are named mol_<N> (vina_redock.py's mol_{i:04d}) and sit
    directly under --results-dir (i.e. --results-dir == vina_redock.py's --out).
  - Folder order (sorted numerically by the trailing integer) corresponds
    1:1, in order, to the row order of the input CSV -- true by construction
    since vina_redock.py names folders from enumerate(reader) over that
    same CSV.
  - Each folder contains docked.pdbqt with Vina's
    "REMARK VINA RESULT:" lines; the first one is the best (lowest/most
    negative) pose score, per Vina's own output ordering -- the same line
    vina_redock.py itself reads for its printed score.
"""

import argparse
import csv
import re
import shutil
import sys
from pathlib import Path

MOL_DIR_RE = re.compile(r"^mol_(\d+)$")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--csv", required=True, type=Path, help="Input dock_remaining.csv with SMILES column")
    p.add_argument("--results-dir", default=Path("."), type=Path, help="Dir containing mol_XXXX folders")
    p.add_argument("--threshold", required=True, type=float,
                   help="Docking score cutoff in kcal/mol. Folders with score <= threshold are kept "
                        "(more negative = better binding, e.g. -11.0)")
    p.add_argument("--out-dir", default=Path("top_candidates"), type=Path, help="Destination dir for top candidates")
    p.add_argument("--out-csv", default=Path("redock_out_sorted.csv"), type=Path, help="Output ranked CSV path")
    p.add_argument("--smiles-col", default="SMILES", help="Name of the SMILES column in the input CSV")
    p.add_argument("--copy", action="store_true", help="Copy folders instead of moving them")
    p.add_argument("--dry-run", action="store_true", help="Report what would happen without moving/copying files")
    return p.parse_args()


def load_smiles(csv_path: Path, smiles_col: str) -> list[str]:
    with csv_path.open(newline="") as f:
        reader = csv.DictReader(f)
        if smiles_col not in reader.fieldnames:
            sys.exit(f"error: column '{smiles_col}' not found in {csv_path}. "
                      f"Available columns: {reader.fieldnames}")
        return [row[smiles_col] for row in reader]


def find_mol_dirs(results_dir: Path) -> list[Path]:
    dirs = []
    for child in results_dir.iterdir():
        if child.is_dir():
            m = MOL_DIR_RE.match(child.name)
            if m:
                dirs.append((int(m.group(1)), child))
    dirs.sort(key=lambda t: t[0])
    return [d for _, d in dirs]


def best_vina_score(docked_pdbqt: Path) -> float | None:
    """Parse docked.pdbqt the same way vina_redock.py's own dock() does:
    first line starting with 'REMARK VINA RESULT', 4th whitespace-separated
    field (the affinity). Vina writes poses best-first, so the first match
    is the best pose."""
    if not docked_pdbqt.exists():
        return None
    with docked_pdbqt.open() as f:
        for line in f:
            if line.startswith("REMARK VINA RESULT"):
                return float(line.split()[3])
    return None


def main():
    args = parse_args()

    smiles_list = load_smiles(args.csv, args.smiles_col)
    mol_dirs = find_mol_dirs(args.results_dir)

    if len(smiles_list) != len(mol_dirs):
        print(f"warning: {len(smiles_list)} SMILES rows vs {len(mol_dirs)} mol_XXXX folders "
              f"— zipping to the shorter of the two, in order.", file=sys.stderr)

    rows = []
    for smiles, folder in zip(smiles_list, mol_dirs):
        docked = folder / "docked.pdbqt"
        score = best_vina_score(docked)
        if score is None:
            print(f"warning: no Vina result found in {docked}, skipping {folder.name}", file=sys.stderr)
            continue
        rows.append({"docking": score, "folder": folder.name, "smiles": smiles, "path": folder})

    # Keep folders that meet/beat the threshold (more negative = better)
    passing = [r for r in rows if r["docking"] <= args.threshold]

    # Rank best (most negative) first
    passing.sort(key=lambda r: r["docking"])

    if args.dry_run:
        print(f"[dry-run] {len(passing)} / {len(rows)} folders pass threshold <= {args.threshold}")
        for r in passing:
            print(f"  {r['docking']:.2f}  {r['folder']}  {r['smiles']}")
        return

    args.out_dir.mkdir(parents=True, exist_ok=True)
    for r in passing:
        dest = args.out_dir / r["folder"]
        if dest.exists():
            print(f"warning: {dest} already exists, skipping move for {r['folder']}", file=sys.stderr)
            continue
        if args.copy:
            shutil.copytree(r["path"], dest)
        else:
            shutil.move(str(r["path"]), str(dest))

    with args.out_csv.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["docking", "folder", "smiles"])
        for r in passing:
            writer.writerow([f"{r['docking']:.2f}", r["folder"], r["smiles"]])

    action = "copied" if args.copy else "moved"
    print(f"{len(passing)} / {len(rows)} folders passed threshold <= {args.threshold}; "
          f"{action} to {args.out_dir}/, ranked CSV written to {args.out_csv}")


if __name__ == "__main__":
    main()
