#!/usr/bin/env python
"""Shuffle and split a .smi file into train/validation sets."""
import argparse, os, random
ap = argparse.ArgumentParser()
ap.add_argument("smi"); ap.add_argument("--val-frac", type=float, default=0.15)
ap.add_argument("--seed", type=int, default=42)
args = ap.parse_args()
lines = [ln for ln in open(args.smi) if ln.strip()]
random.Random(args.seed).shuffle(lines)
n_val = max(1, round(len(lines) * args.val_frac))
val, train = lines[:n_val], lines[n_val:]
base = os.path.splitext(args.smi)[0]
open(f"{base}_train.smi", "w").writelines(train)
open(f"{base}_val.smi", "w").writelines(val)
print(f"total={len(lines)} train={len(train)} val={len(val)}")
