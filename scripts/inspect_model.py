#!/usr/bin/env python
"""Dump non-tensor metadata from a REINVENT .model checkpoint.

Usage: python inspect_model.py models/some.model
Shows top-level keys and any small scalar/string/dict values (skips weight
tensors), which is where REINVENT stashes version, model_type, network params,
vocabulary, and sometimes training info.
"""
import sys, torch

ckpt = torch.load(sys.argv[1], map_location="cpu", weights_only=False)

def show(obj, prefix="", depth=0):
    if depth > 3:
        return
    if isinstance(obj, dict):
        for k, v in obj.items():
            kp = f"{prefix}.{k}" if prefix else str(k)
            if hasattr(v, "shape"):                      # a tensor -> just note shape
                print(f"{kp}: <tensor {tuple(v.shape)}>")
            elif isinstance(v, dict):
                # if it's a state_dict (all tensors), summarize instead of recursing
                if v and all(hasattr(x, "shape") for x in v.values()):
                    print(f"{kp}: <state_dict, {len(v)} tensors>")
                else:
                    print(f"{kp}/")
                    show(v, kp, depth + 1)
            elif isinstance(v, (str, int, float, bool, type(None))):
                print(f"{kp}: {v!r}")
            elif isinstance(v, (list, tuple)):
                preview = str(v)[:200]
                print(f"{kp}: {type(v).__name__}[{len(v)}] {preview}")
            else:
                print(f"{kp}: {type(v).__name__}")
    else:
        print(f"{prefix}: {type(obj).__name__}")

print(f"=== top-level type: {type(ckpt).__name__} ===")
show(ckpt)
