#!/usr/bin/env python
"""Convert SDF file(s) to a standardized SMILES file for REINVENT4 mol2mol TL."""
import argparse, glob, os, sys
from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize
RDLogger.DisableLog("rdApp.*")

def standardize(mol, keep_stereo=True):
    if mol is None:
        return None
    try:
        mol = rdMolStandardize.Cleanup(mol)            # sanitize, disconnect metals, normalize
        mol = rdMolStandardize.FragmentParent(mol)     # keep largest fragment (strips salts)
        mol = rdMolStandardize.Uncharger().uncharge(mol)
        if mol is None or mol.GetNumAtoms() == 0:
            return None
        return Chem.MolToSmiles(mol, isomericSmiles=keep_stereo, canonical=True)
    except Exception:
        return None

def iter_sdf_paths(inputs):
    for p in inputs:
        if os.path.isdir(p):
            yield from sorted(glob.glob(os.path.join(p, "*.sdf")))
        else:
            yield p

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inputs", nargs="+", help="SDF file(s) or a directory of SDFs")
    ap.add_argument("-o", "--output", default="compounds.smi")
    ap.add_argument("--no-stereo", action="store_true")
    args = ap.parse_args()

    seen, n_read, n_fail = {}, 0, 0
    for path in iter_sdf_paths(args.inputs):
        for mol in Chem.SDMolSupplier(path, removeHs=False, sanitize=True):
            n_read += 1
            name = mol.GetProp("_Name").strip() if (mol and mol.HasProp("_Name")) else f"mol{n_read}"
            smi = standardize(mol, keep_stereo=not args.no_stereo)
            if smi is None:
                n_fail += 1
                print(f"  [skip] could not standardize: {name}", file=sys.stderr)
                continue
            seen.setdefault(smi, name)
    with open(args.output, "w") as fh:
        for smi, name in seen.items():
            fh.write(f"{smi}\t{name}\n")
    print(f"read={n_read} failed={n_fail} unique_written={len(seen)} -> {args.output}")

if __name__ == "__main__":
    main()
