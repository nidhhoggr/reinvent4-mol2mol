#!/usr/bin/env python3
"""
vina_redock.py -- deterministic standalone AutoDock Vina docking that mirrors
the exact pipeline DockStream used inside your REINVENT4 run.

Every stochastic stage is pinned to seed 42, identical to DockStream:
  1. SMILES -> 3D   : RDKit EmbedMolecule(randomSeed=42, useRandomCoords=True)
                      -> UFFOptimizeMolecule(maxIters=600) -> AddHs(addCoords=True)
                      -> write PDB  (DockStream writes PDB, NOT sdf, then converts)
  2. PDB -> PDBQT   : obabel <pdb> -opdbqt -O<out> --partialcharge gasteiger
                      (no -xr, so the ligand keeps its rotatable-bond tree)
  3. dock           : vina ... --cpu 1 --seed 42 --exhaustiveness N

Run it with a Python that has RDKit (e.g. the DockStream or reinvent4 conda env),
and point --obabel / --vina at the SAME binaries the container used, or the
scores will not match your run (version differences shift atom typing / search).

Parallelism note:
  Each individual docking stays pinned to --cpu 1 / --seed 42 so it is
  bit-for-bit reproducible against your DockStream run. To use multiple cores,
  parallelize ACROSS ligands in the --csv path with --workers N (one ligand per
  core). Do NOT raise --cpu instead: Vina's multithreaded search is not
  reproducible across thread counts even with a fixed seed.

Single molecule:
  python vina_redock.py "CCO..." \
    --receptor /workspace/docking_setup/receptor.pdbqt \
    --center 1.56 2.56 18.4 --size 22 22 22 --exhaustiveness 8 \
    --vina /opt/conda/envs/reinvent4/bin/vina \
    --obabel /opt/conda/envs/DockStream/bin/obabel \
    --out /workspace/results/redock

Batch from a CSV (uses the SMILES column), 10 ligands at a time:
  python vina_redock.py --csv /workspace/results/combined.csv --smiles-col SMILES \
    --workers 10 ...
"""
import argparse
import csv
import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed

from rdkit import Chem
from rdkit.Chem import AllChem

SEED = 42  # DockStream pins this everywhere; do not change if you want to match the run


def embed(smiles: str) -> Chem.Mol:
    """SMILES -> optimized 3D RDKit mol, identical to DockStream's RDkit pool."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles}")
    code = AllChem.EmbedMolecule(mol, randomSeed=SEED, useRandomCoords=True)
    if code == -1:
        raise RuntimeError("RDKit embedding failed (returned -1) -- "
                           "this is the same failure mode as the ~1/128 fails in your run")
    AllChem.UFFOptimizeMolecule(mol, maxIters=600)
    mol = Chem.AddHs(mol, addCoords=True)
    return mol


def to_pdbqt(mol: Chem.Mol, obabel: str, workdir: str) -> str:
    """Optimized 3D mol -> PDB -> PDBQT via OpenBabel + Gasteiger charges."""
    pdb = os.path.join(workdir, "lig.pdb")
    pdbqt = os.path.join(workdir, "lig.pdbqt")
    Chem.MolToPDBFile(mol, pdb)
    subprocess.run(
        [obabel, pdb, "-opdbqt", f"-O{pdbqt}", "--partialcharge", "gasteiger"],
        check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )
    if not os.path.exists(pdbqt) or os.path.getsize(pdbqt) == 0:
        raise RuntimeError("OpenBabel produced no PDBQT")
    return pdbqt


def dock(pdbqt: str, receptor: str, center, size, vina: str,
         out_pdbqt: str, exhaustiveness: int, cpu: int) -> float:
    """Run Vina deterministically and return the best (mode 1) affinity.

    --cpu is intentionally 1: this keeps the search single-threaded and therefore
    reproducible against the DockStream run. Parallelism comes from running many
    of these at once (see --workers), not from threading a single dock.
    """
    cx, cy, cz = center
    sx, sy, sz = size
    cmd = [
        vina,
        "--receptor", receptor,
        "--ligand", pdbqt,
        "--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
        "--size_x", str(sx), "--size_y", str(sy), "--size_z", str(sz),
        "--cpu", str(cpu),
        "--seed", str(SEED),
        "--exhaustiveness", str(exhaustiveness),
        "--out", out_pdbqt,
    ]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    with open(out_pdbqt) as fh:
        for line in fh:
            if line.startswith("REMARK VINA RESULT"):
                return float(line.split()[3])
    raise RuntimeError("No 'REMARK VINA RESULT' found in Vina output")


def run_one(smiles: str, name: str, args, out_dir: str):
    work = os.path.join(out_dir, name)
    os.makedirs(work, exist_ok=True)
    out_pdbqt = os.path.join(work, "docked.pdbqt")
    mol = embed(smiles)
    pdbqt = to_pdbqt(mol, args.obabel, work)
    score = dock(pdbqt, args.receptor, args.center, args.size,
                 args.vina, out_pdbqt, args.exhaustiveness, args.cpu)
    return score, out_pdbqt


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("smiles", nargs="?", help="A single SMILES string to dock")
    p.add_argument("--csv", help="CSV file to batch-dock instead of a single SMILES")
    p.add_argument("--smiles-col", default="SMILES", help="SMILES column name in --csv")
    p.add_argument("--receptor", required=True, help="receptor .pdbqt path")
    p.add_argument("--center", nargs=3, type=float, required=True,
                   metavar=("X", "Y", "Z"), help="search box center")
    p.add_argument("--size", nargs=3, type=float, default=[22, 22, 22],
                   metavar=("X", "Y", "Z"), help="search box size (default 22 22 22)")
    p.add_argument("--exhaustiveness", type=int, default=8,
                   help="MATCH your dockstream_config.json value (Vina default 8)")
    p.add_argument("--cpu", type=int, default=1,
                   help="tell vina to use more cpus. don't use for parellel docking(default 1)")
    p.add_argument("--workers", type=int, default=1,
                   help="parallel dockings at once for --csv (each still uses --cpu 1). "
                        "Set to your core count, e.g. 10. No effect on single-SMILES mode.")
    p.add_argument("--vina", default="vina", help="path to AutoDock Vina 1.1.2 binary")
    p.add_argument("--obabel", default="obabel", help="path to obabel binary")
    p.add_argument("--out", default="./redock", help="output directory")
    args = p.parse_args()

    os.makedirs(args.out, exist_ok=True)

    if args.csv:
        with open(args.csv) as fh:
            reader = csv.DictReader(fh)
            if args.smiles_col not in reader.fieldnames:
                sys.exit(f"Column '{args.smiles_col}' not in CSV. "
                         f"Found: {reader.fieldnames}")
            jobs = [(row[args.smiles_col], f"mol_{i:04d}")
                    for i, row in enumerate(reader)]

        results = []
        with ProcessPoolExecutor(max_workers=max(1, args.workers)) as pool:
            futs = {pool.submit(run_one, smi, name, args, args.out): (smi, name)
                    for smi, name in jobs}
            for fut in as_completed(futs):
                smi, name = futs[fut]
                try:
                    score, pose = fut.result()
                    results.append((name, score, smi))
                    print(f"{score:8.2f}  {name}  {smi}", flush=True)
                except Exception as e:
                    print(f"   FAIL  {name}  {smi}  -> {e}", flush=True)

        # Re-emit a stable, input-ordered summary (as_completed prints in finish
        # order). Comment this block out if you don't want the second pass.
        if results:
            print("\n# ---- sorted by molecule index ----", flush=True)
            for name, score, smi in sorted(results, key=lambda r: r[0]):
                print(f"{score:8.2f}  {name}  {smi}", flush=True)
    else:
        if not args.smiles:
            sys.exit("Provide a SMILES string or --csv")
        score, pose = run_one(args.smiles, "single", args, args.out)
        print(f"\nVina affinity: {score:.2f} kcal/mol")
        print(f"Docked pose:   {pose}")


if __name__ == "__main__":
    main()
