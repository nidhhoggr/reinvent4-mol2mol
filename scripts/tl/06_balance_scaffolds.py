#!/usr/bin/env python
"""Rebalance a .smi set by capping over-represented chemotypes.

Prevents scaffold collapse in TL: when one chemotype dominates the training set
(e.g. a docking-favored series merged in over iterations), the model overfits to
it. This thins any group above a fraction ceiling, keeping a MaxMin-diverse
representative subset, and leaves every other group untouched. It only ever
REMOVES molecules, never adds.

The grouping scheme (--group-by) decides what counts as "one chemotype":
  scaffold  Bemis-Murcko scaffold (exact ring systems + linkers). Splits
            isomeric variants (quinoline vs quinazoline) into separate groups.
  generic   Bemis-Murcko framework (atoms->C, bonds->single). Merges some
            variants but still splits different ring topologies.
  butina    Fingerprint similarity clusters (ECFP4 + Butina). Groups a whole
            similar family into one cluster regardless of ring isomerism --
            the right choice when one chemotype hides across several scaffolds.

Usage:
    # inspect with similarity clustering (default), write nothing:
    python 06_balance_scaffolds.py configs/compounds.smi --report-only

    # compare groupings on the same set:
    python 06_balance_scaffolds.py configs/compounds.smi --report-only --group-by scaffold
    python 06_balance_scaffolds.py configs/compounds.smi --report-only --group-by butina --cluster-cutoff 0.5

    # write a rebalanced set, 20% ceiling, protecting references:
    python 06_balance_scaffolds.py configs/compounds.smi \
        -o configs/compounds_balanced.smi --group-by butina --max-frac 0.2 \
        --protect configs/references.smi
"""
import argparse
import sys
from collections import defaultdict

from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.SimDivFilters.rdSimDivPickers import MaxMinPicker
from rdkit.ML.Cluster import Butina

RDLogger.DisableLog("rdApp.*")


def canon(smi):
    m = Chem.MolFromSmiles(smi)
    return Chem.MolToSmiles(m) if m else None


def scaffold_of(smi, generic=False):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    try:
        scaf = MurckoScaffold.GetScaffoldForMol(m)
        if scaf is None or scaf.GetNumAtoms() == 0:
            return "(acyclic)"
        if generic:
            scaf = MurckoScaffold.MakeScaffoldGeneric(scaf)
        return Chem.MolToSmiles(scaf)
    except Exception:
        return None


def read_smi(path):
    out = []
    for ln in open(path):
        s = ln.rstrip("\n")
        if s.strip():
            parts = s.split()
            out.append((parts[0], parts[1] if len(parts) > 1 else "", s))
    return out


def load_protect(path):
    prot = set()
    if path:
        for smi, _, _ in read_smi(path):
            c = canon(smi)
            if c:
                prot.add(c)
    return prot


def _fp(smi):
    m = Chem.MolFromSmiles(smi)
    return AllChem.GetMorganFingerprintAsBitVect(m, 2, 2048) if m else None


def maxmin_pick(smiles_list, k, seed=42):
    fps = [_fp(s) for s in smiles_list]
    valid = [i for i, f in enumerate(fps) if f is not None]
    if len(valid) <= k:
        return valid
    vfps = [fps[i] for i in valid]

    def dist(i, j):
        return 1.0 - DataStructs.TanimotoSimilarity(vfps[i], vfps[j])

    picked = MaxMinPicker().LazyPick(dist, len(vfps), k, [], seed)
    return [valid[p] for p in picked]


def butina_clusters(rows, sim_cutoff):
    """Cluster rows by ECFP4 Tanimoto. Returns (clusters, n_bad) where each
    cluster is a list of row indices, centroid first, largest cluster first."""
    fps, orig = [], []
    n_bad = 0
    for i, (smi, _, _) in enumerate(rows):
        f = _fp(smi)
        if f is None:
            n_bad += 1
        else:
            fps.append(f)
            orig.append(i)
    n = len(fps)
    dists = []
    for i in range(n):
        if i:
            sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
            dists.extend(1.0 - s for s in sims)
    clusters = Butina.ClusterData(dists, n, 1.0 - sim_cutoff, isDistData=True)
    return [[orig[j] for j in cl] for cl in clusters], n_bad


def compute_groups(rows, scheme, sim_cutoff):
    """Return (groups: label->[row idx], labels: label->display, n_bad)."""
    groups, labels = {}, {}
    if scheme in ("scaffold", "generic"):
        tmp = defaultdict(list)
        n_bad = 0
        for i, (smi, _, _) in enumerate(rows):
            sc = scaffold_of(smi, generic=(scheme == "generic"))
            if sc is None:
                n_bad += 1
            else:
                tmp[sc].append(i)
        for sc, idxs in tmp.items():
            groups[sc] = idxs
            labels[sc] = sc
        return groups, labels, n_bad
    # butina
    clusters, n_bad = butina_clusters(rows, sim_cutoff)
    for k, cl in enumerate(clusters):
        label = f"cluster_{k}"
        groups[label] = cl
        centroid_scaf = scaffold_of(rows[cl[0]][0]) or "?"
        disp = (centroid_scaf[:40] + "…") if len(centroid_scaf) > 40 else centroid_scaf
        labels[label] = f"cluster {k} · ≈{disp}"
    return groups, labels, n_bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input")
    ap.add_argument("-o", "--output")
    ap.add_argument("--group-by", choices=["scaffold", "generic", "butina"],
                    default="butina", help="how to define a chemotype (default butina)")
    ap.add_argument("--cluster-cutoff", type=float, default=0.5,
                    help="Tanimoto similarity for butina clustering (default 0.5)")
    ap.add_argument("--max-frac", type=float, default=0.2,
                    help="no group may exceed this fraction of the set (default 0.2)")
    ap.add_argument("--max-count", type=int, help="absolute per-group cap (overrides --max-frac)")
    ap.add_argument("--protect", help="SMILES file of molecules to never drop")
    ap.add_argument("--report-only", action="store_true")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    rows = read_smi(args.input)
    total = len(rows)
    if total == 0:
        sys.exit("empty input")

    protect = load_protect(args.protect)
    groups, labels, n_bad = compute_groups(rows, args.group_by, args.cluster_cutoff)
    cap = args.max_count if args.max_count else max(1, int(args.max_frac * total))
    ordered = sorted(groups.items(), key=lambda kv: -len(kv[1]))

    scheme_name = {"scaffold": "Bemis-Murcko", "generic": "generic frameworks",
                   "butina": f"Butina clusters @ Tanimoto {args.cluster_cutoff}"}[args.group_by]
    print(f"=== chemotype histogram ({scheme_name}) ===")
    print(f"{total} molecules, {len(groups)} groups, cap = {cap} ({100*cap/total:.0f}% of set)\n")
    print(f"{'count':>5} {'%':>4}  group")
    print("-" * 70)
    for lab, idxs in ordered[:25]:
        flag = "  <-- OVER CAP" if len(idxs) > cap else ""
        print(f"{len(idxs):>5} {100*len(idxs)/total:>3.0f}%  {labels[lab]}{flag}")
    if len(ordered) > 25:
        tail = sum(len(i) for _, i in ordered[25:])
        print(f"  ... +{len(ordered)-25} more groups ({tail} molecules)")

    over = [(l, i) for l, i in ordered if len(i) > cap]
    print(f"\n{len(over)} group(s) over cap.")
    if n_bad:
        print(f"({n_bad} molecules unparseable / skipped)")

    if args.report_only:
        return
    if not args.output:
        sys.exit("\nProvide -o OUTPUT to write the rebalanced set (or use --report-only).")

    keep = set()
    for lab, idxs in groups.items():
        if len(idxs) <= cap:
            keep.update(idxs)
            continue
        prot_idx = [i for i in idxs if canon(rows[i][0]) in protect]
        rest = [i for i in idxs if i not in set(prot_idx)]
        keep.update(prot_idx)
        slots = cap - len(prot_idx)
        if slots > 0 and rest:
            rest_smi = [rows[i][0] for i in rest]
            for p in maxmin_pick(rest_smi, min(slots, len(rest)), seed=args.seed):
                keep.add(rest[p])

    kept_rows = [rows[i] for i in sorted(keep)]
    with open(args.output, "w") as fh:
        for _, _, raw in kept_rows:
            fh.write(raw + "\n")

    # report largest chemotype in the KEPT set, regrouped the same way
    fgroups, _, _ = compute_groups(kept_rows, args.group_by, args.cluster_cutoff)
    top_n = max((len(v) for v in fgroups.values()), default=0)
    print(f"\nkept {len(kept_rows)} / {total}  (dropped {total - len(kept_rows)})")
    print(f"largest chemotype now {top_n} ({100*top_n/len(kept_rows):.0f}% of kept set)")
    print(f"-> {args.output}")
    print("Re-run check_pairs.py on the output before training.")


if __name__ == "__main__":
    main()
