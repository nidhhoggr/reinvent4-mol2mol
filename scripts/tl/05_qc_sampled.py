#!/usr/bin/env python
"""QC report for a REINVENT4 mol2mol sampling CSV.

Reports validity, uniqueness, novelty vs the TL training set, the
Tanimoto-to-seed distribution, and flags molecules that are either memorized
training compounds or just the seed handed back unchanged.

Usage:
    python qc_sampled.py results/sampled.csv \
        --train configs/compounds_augmented_train.smi \
        -o results/sampled_qc.csv

If your CSV column names differ from the defaults, pass --smiles-col /
--input-col. If the CSV has no input column, pass --seeds <seed.smi> so
Tanimoto-to-seed can be computed against the nearest seed.
"""
import argparse
import csv
import statistics as stats
import sys

from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")


def canon(smi, iso=True):
    if not isinstance(smi, str) or not smi:
        return None
    m = Chem.MolFromSmiles(smi)
    return Chem.MolToSmiles(m, isomericSmiles=iso, canonical=True) if m else None


def fp(smi):
    m = Chem.MolFromSmiles(smi)
    return AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048) if m else None


def load_smiles_file(path):
    out = set()
    with open(path) as fh:
        for ln in fh:
            ln = ln.strip()
            if ln:
                c = canon(ln.split()[0])
                if c:
                    out.add(c)
    return out


def detect_col(header, candidates):
    low = {h.lower(): h for h in header}
    for c in candidates:
        if c.lower() in low:
            return low[c.lower()]
    return None


def pct(n, d):
    return f"{100.0 * n / d:.1f}%" if d else "n/a"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv", help="REINVENT sampling output CSV")
    ap.add_argument("--train", help="TL training SMILES file (for novelty/memorization)")
    ap.add_argument("--seeds", help="seed SMILES file (only if CSV has no input column)")
    ap.add_argument("--smiles-col", help="override generated-SMILES column name")
    ap.add_argument("--input-col", help="override input/seed column name")
    ap.add_argument("-o", "--output", help="write an annotated per-molecule CSV here")
    args = ap.parse_args()

    with open(args.csv, newline="") as fh:
        reader = csv.DictReader(fh)
        header = reader.fieldnames or []
        rows = list(reader)
    if not rows:
        sys.exit("No rows in CSV.")

    smi_col = args.smiles_col or detect_col(header, ["SMILES", "smiles", "Canonical_SMILES"])
    in_col = args.input_col or detect_col(header, ["Input_SMILES", "Input", "input_smiles", "Seed"])
    nll_col = detect_col(header, ["NLL", "nll"])
    if smi_col is None:
        sys.exit(f"Could not find a SMILES column. Columns present: {header}\n"
                 f"Pass --smiles-col explicitly.")

    train = load_smiles_file(args.train) if args.train else set()
    seed_set = load_smiles_file(args.seeds) if args.seeds else set()
    seed_fps = [f for f in (fp(s) for s in seed_set) if f is not None]

    total = len(rows)
    n_valid = 0
    seen = set()
    annotated = []
    tani = []

    for r in rows:
        c = canon(r.get(smi_col, ""))
        valid = c is not None
        if valid:
            n_valid += 1
        is_dup = valid and c in seen
        if valid and not is_dup:
            seen.add(c)

        in_smi = canon(r.get(in_col, "")) if in_col else None

        t = None
        if valid:
            gfp = fp(c)
            if gfp is not None:
                if in_smi:
                    ifp = fp(in_smi)
                    if ifp is not None:
                        t = DataStructs.TanimotoSimilarity(gfp, ifp)
                elif seed_fps:
                    t = max(DataStructs.BulkTanimotoSimilarity(gfp, seed_fps))
        if t is not None:
            tani.append(t)

        memorized = valid and c in train
        equals_seed = valid and ((in_smi is not None and c == in_smi) or (c in seed_set))
        novel = valid and not memorized and not equals_seed

        annotated.append({
            "canonical_smiles": c or "",
            "input_smiles": in_smi or (r.get(in_col, "") if in_col else ""),
            "tanimoto_to_seed": f"{t:.3f}" if t is not None else "",
            "nll": r.get(nll_col, "") if nll_col else "",
            "valid": int(valid),
            "duplicate": int(is_dup),
            "memorized_training": int(memorized),
            "equals_seed": int(equals_seed),
            "novel": int(novel),
        })

    unique = len(seen)
    n_memorized = sum(a["memorized_training"] for a in annotated)
    n_seedcopy = sum(a["equals_seed"] for a in annotated)
    novel_unique = {a["canonical_smiles"] for a in annotated if a["novel"]}
    n_novel_unique = len(novel_unique)

    print(f"=== QC: {args.csv} ===")
    print(f"rows (generated):        {total}")
    print(f"valid:                   {n_valid}  ({pct(n_valid, total)})")
    print(f"unique (valid):          {unique}  ({pct(unique, n_valid)} of valid)")
    if args.train:
        print(f"memorized (in train):    {n_memorized}  ({pct(n_memorized, n_valid)} of valid)")
    else:
        print("memorized (in train):    [skipped — pass --train to enable]")
    print(f"seed returned unchanged: {n_seedcopy}  ({pct(n_seedcopy, n_valid)} of valid)")
    print(f"novel & unique:          {n_novel_unique}  <- your usable pool")

    if tani:
        ts = sorted(tani)
        def q(p):
            return ts[min(len(ts) - 1, int(p * len(ts)))]
        print("\nTanimoto-to-seed (ECFP4):")
        print(f"  min {min(tani):.2f} | 25% {q(0.25):.2f} | median {stats.median(tani):.2f} "
              f"| 75% {q(0.75):.2f} | max {max(tani):.2f} | mean {stats.fmean(tani):.2f}")
        bands = [("identical (~1.0)", lambda x: x >= 0.999),
                 ("very close .8-1", lambda x: 0.8 <= x < 0.999),
                 ("close .6-.8", lambda x: 0.6 <= x < 0.8),
                 ("moderate .4-.6", lambda x: 0.4 <= x < 0.6),
                 ("distant <.4", lambda x: x < 0.4)]
        print("  distribution:")
        for label, fn in bands:
            cnt = sum(1 for x in tani if fn(x))
            print(f"    {label:17s}: {cnt:5d}  ({pct(cnt, len(tani))})")
    else:
        print("\nTanimoto-to-seed: [no input/seed available — pass --seeds or check --input-col]")

    if args.output:
        cols = ["canonical_smiles", "input_smiles", "tanimoto_to_seed", "nll",
                "valid", "duplicate", "memorized_training", "equals_seed", "novel"]
        with open(args.output, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            w.writeheader()
            w.writerows(annotated)
        print(f"\nannotated per-molecule CSV -> {args.output}")


if __name__ == "__main__":
    main()
