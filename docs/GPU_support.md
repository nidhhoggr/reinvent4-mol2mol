## Docker

```yaml
# docker-compose.gpu.yml
services:
  builder:
    build:
      args:
        TORCH_INDEX: ${TORCH_INDEX:-https://download.pytorch.org/whl/cu126}

  tl:    &gpu { deploy: { resources: { reservations: { devices: [{ driver: nvidia, count: all, capabilities: [gpu] }] } } } }
  sample: *gpu
  rl:     *gpu
```

One real gotcha with putting build args **only** in the override file: in CPU mode that file isn't loaded, so `builder.build.args` is absent entirely — which is fine, because your **Dockerfile** still has `ARG TORCH_INDEX=…/cpu` as the fallback. So the resolution chain is:

- **CPU run** (`docker compose up`): GPU file not loaded → no build arg passed → Dockerfile `ARG` default (`/cpu`) wins. ✓
- **GPU run** (`-f … -f docker-compose.gpu.yml`): override supplies `cu126` → Dockerfile `ARG` overridden. ✓

The thing to *not* do is split it: TORCH_INDEX default in the base `builder` and the GPU override in the gpu file. That works but now two files describe the same arg and you have to reason about precedence. Pick one home for it.

``bash
# CPU — gpu file not loaded, Dockerfile ARG default (/cpu) wins
docker compose up tl

# GPU — override passes cu126 at build AND reserves the device at run
docker compose -f docker-compose.yml -f docker-compose.gpu.yml build tl
docker compose -f docker-compose.yml -f docker-compose.gpu.yml up tl

# The env var feeds the ${TORCH_INDEX} substitution in the YAML, overriding the cu126 default. No --build-arg needed at all. And because it's the same TORCH_INDEX variable, you can drop it in your root .env and forget about it.
TORCH_INDEX=https://download.pytorch.org/whl/cu130 \
  docker compose -f docker-compose.yml -f docker-compose.gpu.yml build
```

Compose deep-merges the two, so the GPU file only contributes the device reservation. The base file is byte-for-byte identical between the two modes.

If typing `-f … -f …` every time annoys you and you want it **env-driven** (closer to what you're after), use the special `.env` file's `COMPOSE_FILE` variable:

```
# .env  (project root)
COMPOSE_FILE=docker-compose.yml:docker-compose.gpu.yml
TORCH_INDEX=https://download.pytorch.org/whl/cu130
```
Now a plain `docker compose up` picks up both. Flip to CPU by commenting that line out or setting `COMPOSE_FILE=docker-compose.yml`.

The one catch to keep straight: that `.env` is the special file Compose auto-reads for `${VAR}` substitution and settings like `COMPOSE_FILE` — it is **not** your `host.env`. `host.env` stays doing its job (runtime `OMP_NUM_THREADS` etc. inside the container); `.env` drives Compose's own behavior. They're two different mechanisms that happen to both involve env-style files.

This also pairs cleanly with the build-arg side from before — you can put `TORCH_INDEX` in that same `.env` so one file selects both the CUDA wheel at build and the GPU reservation at run.

## TOML

Finally, set `device = "cuda:0"` in each config (`tl_mol2mol.toml`, `sample_mol2mol.toml`, `staged_learning.toml`). With a CPU build this line errors or silently falls back, so it's the last switch to flip.

---

One unrelated flag while you're in there: `torch==2.12.0` is a very high pin — worth confirming it actually resolves on `cu126` before you commit to a long GPU build, since a nonexistent version will fail Layer 3 regardless of the index. If you want, I can check what the current stable torch is and which CUDA wheel indexes it ships for.

Two honest caveats so the new machine doesn't disappoint you:

The GPU only accelerates the **generator** (the LibInvent RNN). Your **docking is CPU-bound** — DockStream/Vina don't use the GPU — and on this pipeline docking dominates wall-clock time. So a faster GPU speeds sampling and the policy update, but throughput will still be gated by how fast 128 ligands dock per step. A machine with more/faster CPU cores helps there too; bump `number_cores` in `dockstream_config.json` (and you can raise `batch_size` back up once the generator is on GPU).

And match `TORCH_INDEX` to whatever torch version REINVENT pins at the time you rebuild — it's `torch==2.12.0` now, and the cu126 index carries that build, but if the pin moves you'd point at the corresponding CUDA index (cu128, etc.). The Layer 5 line now prints the active torch and CUDA build at build time, so you'll see immediately whether you got the CPU or CUDA wheel.
