#!/usr/bin/env python3
"""Dump REINVENT4 TensorBoard scalars to a wide CSV (step x tag) + print a summary.

Usage:
    python dump_tb.py [LOGDIR] [OUTDIR]

LOGDIR defaults to "tb_rl", OUTDIR to ".". Handles nested/resumed run dirs.
Only dependency is `tensorboard` (already present in the reinvent env).
"""
import csv
import os
import sys

try:
    from tensorboard.backend.event_processing.event_accumulator import EventAccumulator
except ImportError:
    sys.exit("tensorboard not importable -- run inside the reinvent container, "
             "or `pip install tensorboard`.")

# Column order for the wide CSV / summary; unknown tags get appended alphabetically.
PREFERRED = [
    "Average total score", "Loss",
    "docking", "docking (raw)",
    "QED", "QED (raw)",
    "MW", "MW (raw)",
    "Halogens", "Halogens (raw)",
    "NLL/agent", "NLL/prior", "NLL/augmented",
    "Fraction of valid SMILES", "Fraction of duplicate SMILES",
    "Number of unique scaffolds",
]


def find_run_dirs(root):
    runs = []
    for d, _, files in os.walk(root):
        if any(f.startswith("events.out.tfevents") for f in files):
            runs.append(d)
    return sorted(runs)


def load_run(run_dir):
    ea = EventAccumulator(run_dir, size_guidance={"scalars": 0})  # 0 = keep all points
    ea.Reload()
    data = {}  # tag -> {step: value}
    for tag in ea.Tags().get("scalars", []):
        data[tag] = {ev.step: ev.value for ev in ea.Scalars(tag)}
    return data


def order_tags(tags):
    ordered = [t for t in PREFERRED if t in tags]
    return ordered + sorted(t for t in tags if t not in ordered)


def write_wide_csv(path, data):
    tags = order_tags(list(data))
    steps = sorted({s for series in data.values() for s in series})
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["step"] + tags)
        for s in steps:
            row = [s]
            for t in tags:
                v = data[t].get(s)
                row.append(f"{v:.6g}" if v is not None else "")
            w.writerow(row)
    return tags, steps


def summarize(data, tags, steps):
    print(f"  steps {steps[0]}..{steps[-1]}  ({len(steps)} points)")
    width = max(len(t) for t in tags)
    hdr = "  " + "tag".ljust(width) + "   first      last        min        max          delta"
    print(hdr)
    for t in tags:
        series = data[t]
        ks = sorted(series)
        first, last = series[ks[0]], series[ks[-1]]
        vals = list(series.values())
        print(f"  {t.ljust(width)}  {first:9.4g}  {last:9.4g}  "
              f"{min(vals):9.4g}  {max(vals):9.4g}  {last - first:+9.4g}")


def main():
    root = sys.argv[1] if len(sys.argv) > 1 else "tb_rl"
    outdir = sys.argv[2] if len(sys.argv) > 2 else "."
    if not os.path.isdir(root):
        sys.exit(f"logdir not found: {root}")
    runs = find_run_dirs(root)
    if not runs:
        sys.exit(f"no events.out.tfevents.* files under {root}")
    os.makedirs(outdir, exist_ok=True)
    multi = len(runs) > 1

    for run_dir in runs:
        data = load_run(run_dir)
        if not data:
            continue
        name = os.path.relpath(run_dir, root).replace(os.sep, "_")
        if name in (".", ""):
            name = os.path.basename(os.path.abspath(root))
        fname = f"tb_stats_{name}.csv" if multi else "tb_stats.csv"
        out = os.path.join(outdir, fname)
        tags, steps = write_wide_csv(out, data)
        print(f"\n[{run_dir}] -> {out}")
        summarize(data, tags, steps)

    print("\nUpload the tb_stats*.csv here and we can walk the trends.")


if __name__ == "__main__":
    main()
