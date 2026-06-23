#!/usr/bin/env python
"""Compare two sampling sets: overlap (duplicates), what's unique to each, and
scaffold-level overlap. A high overlap fraction across iterations is a direct
sign the model is converging/collapsing rather than exploring.

Usage:
    python compare_samples.py iter1.csv iter2.csv
    python compare_samples.py iter1.smi iter2.smi --smiles-col SMILES -o overlap.csv

Accepts .csv (auto-detects a SMILES column) or .smi (first column).
"""
import argparse, csv, os, sys
from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold
RDLogger.DisableLog("rdApp.*")

def canon(s):
    m = Chem.MolFromSmiles(s) if isinstance(s, str) else None
    return Chem.MolToSmiles(m) if m else None

def scaffold(s):
    m = Chem.MolFromSmiles(s)
    if m is None: return None
    sc = MurckoScaffold.GetScaffoldForMol(m)
    return Chem.MolToSmiles(sc) if sc and sc.GetNumAtoms() else "(acyclic)"

def detect_col(header, cands):
    low = {h.lower(): h for h in header}
    for c in cands:
        if c.lower() in low: return low[c.lower()]
    return None

def load(path, smiles_col=None):
    smis = []
    if path.lower().endswith((".smi", ".txt")):
        for ln in open(path):
            if ln.strip(): smis.append(ln.split()[0])
    else:
        with open(path, newline="") as fh:
            r = csv.DictReader(fh); rows = list(r)
            col = smiles_col or detect_col(r.fieldnames or [], ["SMILES","smiles","Canonical_SMILES"])
            if not col: sys.exit(f"No SMILES column in {path}; columns={r.fieldnames}")
            smis = [row.get(col,"") for row in rows]
    out = set()
    for s in smis:
        c = canon(s)
        if c: out.add(c)
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("file_a"); ap.add_argument("file_b")
    ap.add_argument("--smiles-col")
    ap.add_argument("-o","--output", help="write the shared (duplicate) molecules here")
    args = ap.parse_args()

    A = load(args.file_a, args.smiles_col)
    B = load(args.file_b, args.smiles_col)
    shared = A & B
    only_a, only_b = A - B, B - A
    union = A | B
    na, nb = os.path.basename(args.file_a), os.path.basename(args.file_b)

    print(f"=== sample overlap ===")
    print(f"{na}: {len(A)} unique molecules")
    print(f"{nb}: {len(B)} unique molecules")
    print(f"shared (exact duplicates): {len(shared)}")
    print(f"  = {100*len(shared)/len(A):.1f}% of {na}, "
          f"{100*len(shared)/len(B):.1f}% of {nb}, "
          f"{100*len(shared)/len(union):.1f}% Jaccard")
    print(f"unique to {na}: {len(only_a)}")
    print(f"unique to {nb}: {len(only_b)}")

    # scaffold-level overlap (broader than exact-molecule)
    sa = {scaffold(s) for s in A} - {None}
    sb = {scaffold(s) for s in B} - {None}
    sshared = sa & sb
    print(f"\nscaffold-level: {na} {len(sa)} scaffolds, {nb} {len(sb)}, "
          f"shared {len(sshared)} ({100*len(sshared)/len(sa|sb):.1f}% Jaccard)")
    print("(scaffold overlap >> molecule overlap means same chemotypes, new decorations)")

    if args.output:
        with open(args.output,"w",newline="") as fh:
            w = csv.writer(fh); w.writerow(["shared_smiles"])
            for s in sorted(shared): w.writerow([s])
        print(f"\nshared molecules -> {args.output}")

if __name__ == "__main__":
    main()
