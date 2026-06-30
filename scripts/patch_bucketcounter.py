#!/usr/bin/env python3
"""Patch REINVENT's BucketCounter so its diversity-filter state survives
torch.save/torch.load (i.e. checkpoint resume).

Why: BucketCounter subclasses collections.Counter with a custom
__init__(self, max_size, ...). Counter's default __reduce__ reconstructs via
BucketCounter(dict_of_counts) -> the counts dict gets bound to `max_size` and
the counts are dropped. Result: every checkpoint resume reloads an EMPTY,
max_size-corrupted scaffold memory, silently disabling the anti-collapse filter.

This inserts a correct __reduce__ that rebuilds via (max_size,) and replays the
items, preserving both counts and max_size. Idempotent; fails loudly if the
REINVENT layout has shifted (so it can't silently no-op).

Run at image-build time (editable install -> edit persists). Optional arg = path
to bucket_counter.py; otherwise it searches /opt/reinvent4.
"""
import pathlib
import sys

def find_target(argv):
    if len(argv) > 1:
        return pathlib.Path(argv[1])
    direct = pathlib.Path("/opt/reinvent4/reinvent/runmodes/RL/memories/bucket_counter.py")
    if direct.exists():
        return direct
    hits = list(pathlib.Path("/opt/reinvent4").rglob("memories/bucket_counter.py"))
    if not hits:
        sys.exit("ERROR: bucket_counter.py not found under /opt/reinvent4 "
                 "(REINVENT moved? pass the path explicitly).")
    return hits[0]

def main():
    p = find_target(sys.argv)
    src = p.read_text()

    if "__reduce__" in src:
        print(f"BucketCounter already patched: {p}")
        return

    anchor = "        self.max_size = max_size\n"
    if anchor not in src:
        sys.exit("ERROR: anchor 'self.max_size = max_size' not found in "
                 f"{p} -- REINVENT version changed; update this patch.")

    method = (
        "        self.max_size = max_size\n"
        "\n"
        "    def __reduce__(self):\n"
        "        # PATCH: Counter's default __reduce__ calls BucketCounter(counts),\n"
        "        # binding the counts dict to max_size and dropping the counts, which\n"
        "        # corrupts the diversity filter on every checkpoint resume. Rebuild\n"
        "        # via (max_size,) and replay items so torch.save/load preserves both.\n"
        "        return (self.__class__, (self.max_size,), None, None, iter(self.items()))\n"
    )
    p.write_text(src.replace(anchor, method, 1))
    print(f"Patched BucketCounter.__reduce__ in {p}")

if __name__ == "__main__":
    main()
