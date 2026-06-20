#!/usr/bin/env python
"""Regenerate libinvent_rl.template.toml from the validated libinvent_rl.toml.

Run this every time you change the scoring function (or anything else) in
libinvent_rl.toml, BEFORE launching run_resilient.sh, so the wrapper runs your
real config. It only swaps the few fields the wrapper controls and leaves your
scoring components untouched. Commented (#) lines are never modified.

  python make_template.py /workspace/configs/libinvent_rl.toml \
                          /workspace/configs/libinvent_rl.template.toml
"""
import re, sys

src, dst = sys.argv[1], sys.argv[2]
ROLLING_CHKPT = "/workspace/results/libinvent.chkpt"

done = {k: False for k in ("agent", "csv", "chkpt", "maxsteps", "maxscore", "minsteps")}
out = []

for ln in open(src).read().splitlines():
    if ln.lstrip().startswith("#"):
        out.append(ln); continue
    indent = ln[:len(ln) - len(ln.lstrip())]
    if not done["agent"] and re.match(r"\s*agent_file\s*=", ln):
        out.append(f'{indent}agent_file = "__AGENT__"')
        out.append(f'{indent}use_checkpoint = __USE__')
        done["agent"] = True; continue
    if not done["csv"] and re.match(r"\s*summary_csv_prefix\s*=", ln):
        out.append(re.sub(r"=.*", '= "__CSVPREFIX__"', ln, count=1)); done["csv"] = True; continue
    if not done["chkpt"] and re.match(r"\s*chkpt_file\s*=", ln):
        out.append(re.sub(r"=.*", f'= "{ROLLING_CHKPT}"', ln, count=1)); done["chkpt"] = True; continue
    if not done["maxsteps"] and re.match(r"\s*max_steps\s*=", ln):
        out.append(re.sub(r"=.*", "= __MAXSTEPS__", ln, count=1)); done["maxsteps"] = True; continue
    if not done["maxscore"] and re.match(r"\s*max_score\s*=", ln):
        out.append(re.sub(r"=.*", "= 1.0", ln, count=1)); done["maxscore"] = True; continue
    if not done["minsteps"] and re.match(r"\s*min_steps\s*=", ln):
        out.append(re.sub(r"=.*", "= 0", ln, count=1)); done["minsteps"] = True; continue
    out.append(ln)

open(dst, "w").write("\n".join(out) + "\n")
missing = [k for k, v in done.items() if not v]
print("Wrote", dst)
if missing:
    print("WARNING: did not find/replace:", missing)
