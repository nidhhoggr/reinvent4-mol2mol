"""Dump REINVENT4 TensorBoard scalars to a wide CSV (step x tag) + print a summary.

Usage:
    python dump_tb.py [LOGDIR] [OUTDIR] [--emit-tb]

LOGDIR defaults to "results", OUTDIR to ".".

Handles the resilient-loop layout: run_resilient.sh writes each chunk's TB to
results/run_<done_steps>_tb, and REINVENT restarts its step counter at 1 every
chunk. This script reads the <done_steps> offset out of each dir name and
reconstructs the GLOBAL step (offset + local step), so chunked runs stitch into
one continuous series instead of N disconnected segments. A single non-chunked
run (no run_<N> in the path) just dumps as-is.

--emit-tb additionally writes OUTDIR/tb_merged/ : a single fresh events file
holding the stitched, globally-renumbered series, so you can view ONE continuous
curve in TensorBoard instead of N overlaid per-chunk runs:
    tensorboard --logdir OUTDIR/tb_merged --host 0.0.0.0 --port 6006

Dependencies: `tensorboard` (read side, always); `torch` (only for --emit-tb).
Both are present in the reinvent env.
"""
import csv
import os
import re
import shutil
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

CHUNK_RE = re.compile(r"run_(\d+)")   # run_resilient.sh encodes the chunk's start step here


def find_run_dirs(root):
    runs = []
    for d, _, files in os.walk(root):
        if any(f.startswith("events.out.tfevents") for f in files):
            runs.append(d)
    return sorted(runs)


def chunk_offset(run_dir):
    """Pull the global-step offset out of a run_<N>_tb dir name; 0 if absent."""
    m = CHUNK_RE.search(os.path.basename(run_dir)) or CHUNK_RE.search(run_dir)
    return int(m.group(1)) if m else 0


def load_run(run_dir):
    ea = EventAccumulator(run_dir, size_guidance={"scalars": 0})  # 0 = keep all points
    ea.Reload()
    data = {}  # tag -> {step: value}
    for tag in ea.Tags().get("scalars", []):
        data[tag] = {ev.step: ev.value for ev in ea.Scalars(tag)}
    return data


def merge_runs(run_dirs):
    """Stitch dirs into one series keyed by global step = offset + local step.
    Dirs are processed in offset order so any boundary-step overlap keeps the
    later chunk's value (last write wins)."""
    merged = {}  # tag -> {global_step: value}
    for rd in sorted(run_dirs, key=chunk_offset):
        off = chunk_offset(rd)
        for tag, series in load_run(rd).items():
            dst = merged.setdefault(tag, {})
            for step, val in series.items():
                dst[off + step] = val
    return merged


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


def emit_tb(data, out_dir):
    """Write the stitched series to a single fresh events file with global steps."""
    try:
        from torch.utils.tensorboard import SummaryWriter
    except ImportError:
        print("  (--emit-tb skipped: torch.utils.tensorboard not importable)")
        return None
    # Start clean so we don't stack a second events file with overlapping steps.
    if os.path.basename(os.path.normpath(out_dir)) == "tb_merged" and os.path.isdir(out_dir):
        shutil.rmtree(out_dir)
    writer = SummaryWriter(out_dir)
    for tag, series in data.items():
        for step, val in sorted(series.items()):
            writer.add_scalar(tag, val, step)
    writer.close()
    return out_dir


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
    argv = sys.argv[1:]
    flags = {a for a in argv if a.startswith("-")}
    pos = [a for a in argv if not a.startswith("-")]
    root = pos[0] if len(pos) > 0 else "results"
    outdir = pos[1] if len(pos) > 1 else "."
    emit = "--emit-tb" in flags

    if not os.path.isdir(root):
        sys.exit(f"logdir not found: {root}")
    runs = find_run_dirs(root)
    if not runs:
        sys.exit(f"no events.out.tfevents.* files under {root}")
    os.makedirs(outdir, exist_ok=True)

    chunked = [r for r in runs if CHUNK_RE.search(os.path.basename(r)) or CHUNK_RE.search(r)]
    if chunked:
        data = merge_runs(chunked)
        label = f"stitched {len(chunked)} chunk(s)"
        skipped = [r for r in runs if r not in chunked]
    else:
        data = load_run(runs[0]) if len(runs) == 1 else merge_runs(runs)
        label = runs[0] if len(runs) == 1 else root
        skipped = []

    out = os.path.join(outdir, "tb_stats.csv")
    tags, steps = write_wide_csv(out, data)
    print(f"\n{label} -> {out}")
    summarize(data, tags, steps)
    for r in skipped:
        print(f"(ignored non-chunk event dir: {r})")

    if emit:
        merged_dir = os.path.join(outdir, "tb_merged")
        res = emit_tb(data, merged_dir)
        if res:
            print(f"\nmerged events -> {res}")
            print(f"  tensorboard --logdir {res} --host 0.0.0.0 --port 6006")

    print("\nUpload tb_stats.csv here and we can walk the trends.")


if __name__ == "__main__":
    main()
