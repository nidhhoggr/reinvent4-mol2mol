#!/usr/bin/env python
"""Report how many (source, target) pairs survive each Tanimoto threshold."""
import sys
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem
RDLogger.DisableLog("rdApp.*")

smis = [ln.split()[0] for ln in open(sys.argv[1]) if ln.strip()]
mols = [m for m in (Chem.MolFromSmiles(s) for s in smis) if m is not None]
fps = [AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048) for m in mols]
n = len(fps)
print(f"{n} valid molecules -> {n*(n-1)} ordered pairs possible (excl. self)\n")
print(f"{'threshold':>9} | {'pairs >= thr':>12} | {'sources w/>=1 target':>20}")
print("-" * 48)
for thr in (0.3, 0.4, 0.5, 0.6, 0.7, 0.8):
    pairs, sources = 0, set()
    for i in range(n):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps)
        for j in range(n):
            if i != j and sims[j] >= thr:
                pairs += 1; sources.add(i)
    print(f"{thr:>9.1f} | {pairs:>12} | {len(sources):>20}")
