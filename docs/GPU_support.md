The good news is that the fix is small and *doesn't* involve wrangling CUDA dependencies, because the CUDA libraries are bundled inside the torch wheel. So a properly built GPU image needs **no CUDA toolkit in the image and no CUDA install on the host** — the host just needs an NVIDIA driver plus the NVIDIA Container Toolkit (standard for any GPU + Docker setup), and you run with `--gpus all`.

This is portable by parameterizing the torch index with a build arg, so the *same* Dockerfile builds CPU or GPU without editing it. Here's that version:So on the new machine, the full GPU path is:

```bash
# build the GPU image (one flag difference)
docker build -t reinvent4-gpu --build-arg TORCH_INDEX=https://download.pytorch.org/whl/cu126 .

# host needs the NVIDIA Container Toolkit once (not a per-project dependency):
#   sudo apt-get install -y nvidia-container-toolkit && sudo nvidia-ctk runtime configure --runtime=docker && sudo systemctl restart docker

# run with the GPU exposed
docker run --gpus all -it -v /path/to/workspace:/workspace reinvent4-gpu
```

Then inside, two checks and one config change:
```bash
nvidia-smi                                                  # GPU visible?
python -c "import torch; print(torch.cuda.is_available())"  # expect True
```
and set `device = "cuda:0"` in `libinvent_rl.toml`.

Why this avoids the dependency-wrangling you're worried about: the cu126 torch wheel ships its own CUDA runtime, so you never install a CUDA toolkit in the image or on the host — the only host-side piece is the driver + container toolkit, which is the standard cost of running any GPU container. And a "more recent" machine sidesteps the exact thing that bit you before: the GTX 970 was Maxwell (sm_52), which modern torch dropped; any current card (RTX 20-series / T4 and newer) is supported, so it'll just work.

Two honest caveats so the new machine doesn't disappoint you:

The GPU only accelerates the **generator** (the LibInvent RNN). Your **docking is CPU-bound** — DockStream/Vina don't use the GPU — and on this pipeline docking dominates wall-clock time. So a faster GPU speeds sampling and the policy update, but throughput will still be gated by how fast 128 ligands dock per step. A machine with more/faster CPU cores helps there too; bump `number_cores` in `dockstream_config.json` (and you can raise `batch_size` back up once the generator is on GPU).

And match `TORCH_INDEX` to whatever torch version REINVENT pins at the time you rebuild — it's `torch==2.12.0` now, and the cu126 index carries that build, but if the pin moves you'd point at the corresponding CUDA index (cu128, etc.). The Layer 5 line now prints the active torch and CUDA build at build time, so you'll see immediately whether you got the CPU or CUDA wheel.
