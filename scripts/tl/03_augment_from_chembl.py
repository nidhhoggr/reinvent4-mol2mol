#!/usr/bin/env python
"""Augment a .smi file with ChEMBL similarity hits to build congeneric clusters
for REINVENT4 mol2mol transfer learning.

For each seed SMILES it queries ChEMBL's similarity endpoint, standardizes the
hits, caps how many it keeps per seed, dedupes against everything, and writes a
new .smi with the originals first followed by the analogs.

Usage:
    python augment_from_chembl.py compounds.smi -o augmented_compounds.smi \
        --similarity 70 --max-per-seed 25

ChEMBL endpoint:
    GET https://www.ebi.ac.uk/chembl/api/data/similarity/{smiles}/{int_similarity}?format=json
    int_similarity is a Tanimoto percentage in [40, 100].
"""
import argparse
import sys
import time
from urllib.parse import quote

import requests
from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize

RDLogger.DisableLog("rdApp.*")

CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"
HEADERS = {"User-Agent": "reinvent4-mol2mol-augment/1.0 (similarity search)"}


# ---------- chemistry helpers ----------

def standardize(smi, keep_stereo=True):
    """Desalt -> neutralize -> canonicalize. Returns canonical SMILES or None."""
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    try:
        mol = rdMolStandardize.Cleanup(mol)
        mol = rdMolStandardize.FragmentParent(mol)
        mol = rdMolStandardize.Uncharger().uncharge(mol)
        if mol is None or mol.GetNumAtoms() == 0:
            return None
        return Chem.MolToSmiles(mol, isomericSmiles=keep_stereo, canonical=True)
    except Exception:
        return None


def heavy_atoms(smi):
    m = Chem.MolFromSmiles(smi)
    return m.GetNumHeavyAtoms() if m else 0


# ---------- ChEMBL fetch (the only networked part) ----------

def chembl_similarity(smi, sim_pct, max_per_seed, session, sleep=0.3, timeout=30):
    """Yield (chembl_id, similarity_float, canonical_smiles) for hits to `smi`.

    Follows pagination until max_per_seed results or no next page. Raises
    requests exceptions only on hard failures; a 404 (no hits) yields nothing.
    """
    url = f"{CHEMBL_BASE}/similarity/{quote(smi, safe='')}/{sim_pct}?format=json&limit=50"
    fetched = 0
    while url and fetched < max_per_seed:
        resp = session.get(url, headers=HEADERS, timeout=timeout)
        if resp.status_code == 404:
            return  # ChEMBL returns 404 when nothing clears the threshold
        resp.raise_for_status()
        data = resp.json()
        for mol in data.get("molecules", []):
            structs = mol.get("molecule_structures") or {}
            cs = structs.get("canonical_smiles")
            if not cs:
                continue
            try:
                sim = float(mol.get("similarity", "0"))
            except (TypeError, ValueError):
                sim = 0.0
            yield mol.get("molecule_chembl_id", "?"), sim, cs
            fetched += 1
            if fetched >= max_per_seed:
                return
        nxt = (data.get("page_meta") or {}).get("next")
        url = f"https://www.ebi.ac.uk{nxt}" if nxt else None
        if url:
            time.sleep(sleep)


# ---------- main ----------

def read_seeds(path):
    seeds = []
    for ln in open(path):
        ln = ln.strip()
        if not ln:
            continue
        parts = ln.split()
        smi = parts[0]
        name = parts[1] if len(parts) > 1 else f"seed{len(seeds)+1}"
        seeds.append((smi, name))
    return seeds


def run(args, fetcher=chembl_similarity):
    seeds = read_seeds(args.input)
    session = requests.Session()

    # ordered: canonical_smiles -> label. Seeds inserted first so they win on dedup.
    out = {}
    for smi, name in seeds:
        std = standardize(smi, keep_stereo=not args.no_stereo)
        if std:
            out.setdefault(std, name)
        else:
            print(f"  [warn] seed failed to standardize: {name}", file=sys.stderr)

    n_seed = len(out)
    per_seed_counts = {}
    empty_seeds = []

    for smi, name in seeds:
        std_seed = standardize(smi, keep_stereo=not args.no_stereo)
        if std_seed is None:
            continue
        added = 0
        try:
            for cid, sim, hit_smi in fetcher(
                std_seed, args.similarity, args.max_per_seed, session, sleep=args.sleep
            ):
                std = standardize(hit_smi, keep_stereo=not args.no_stereo)
                if std is None:
                    continue
                ha = heavy_atoms(std)
                if ha < args.min_heavy or ha > args.max_heavy:
                    continue
                if std in out:          # already have it (seed or prior hit)
                    continue
                out[std] = f"{cid}|from:{name}|sim:{sim/100:.2f}"
                added += 1
        except requests.RequestException as e:
            print(f"  [error] ChEMBL query failed for {name}: {e}", file=sys.stderr)
        per_seed_counts[name] = added
        if added == 0:
            empty_seeds.append(name)
        print(f"  {name}: +{added} analogs", file=sys.stderr)

    with open(args.output, "w") as fh:
        for smi, label in out.items():
            fh.write(f"{smi}\t{label}\n")

    n_total = len(out)
    print(
        f"\nseeds={n_seed} analogs_added={n_total - n_seed} total_unique={n_total} "
        f"-> {args.output}"
    )
    if empty_seeds:
        print(f"seeds with 0 analogs ({len(empty_seeds)}): {', '.join(empty_seeds)}",
              file=sys.stderr)
    print("Now re-run check_pairs.py on the output to see the threshold table fill in.",
          file=sys.stderr)


def parse_args(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("input", help="seed .smi (first column SMILES)")
    ap.add_argument("-o", "--output", default="augmented_compounds.smi")
    ap.add_argument("--similarity", type=int, default=70,
                    help="ChEMBL Tanimoto %% threshold, 40-100 (default 70)")
    ap.add_argument("--max-per-seed", type=int, default=25,
                    help="cap analogs kept per seed (default 25)")
    ap.add_argument("--min-heavy", type=int, default=5)
    ap.add_argument("--max-heavy", type=int, default=70)
    ap.add_argument("--sleep", type=float, default=0.3,
                    help="seconds between paged requests (be polite)")
    ap.add_argument("--no-stereo", action="store_true")
    args = ap.parse_args(argv)
    if not (40 <= args.similarity <= 100):
        ap.error("--similarity must be between 40 and 100")
    return args


if __name__ == "__main__":
    run(parse_args())
