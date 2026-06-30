FROM condaforge/miniforge3:latest

# ============================================================================
# Build for CPU (default) or GPU by overriding TORCH_INDEX at build time:
#
#   CPU:  docker build -t reinvent4 .
#   GPU:  docker build -t reinvent4-gpu \
#             --build-arg TORCH_INDEX=https://download.pytorch.org/whl/cu126 .
#
# The CUDA torch wheel bundles its own CUDA runtime, so the IMAGE needs no CUDA
# toolkit. To actually use the GPU at runtime the HOST needs only:
#   - a recent NVIDIA driver
#   - the NVIDIA Container Toolkit (nvidia-ctk)
# and you launch the container with `--gpus all`. Then set device = "cuda:0"
# in your REINVENT TOML.
# ============================================================================
ARG TORCH_INDEX=https://download.pytorch.org/whl/cpu

# Layer 1: System packages
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    git \
    openbabel \
    && rm -rf /var/lib/apt/lists/*

# Layer 2: Create REINVENT4 conda env
RUN mamba create -n reinvent4 python=3.11 -y

# Layer 3: PyTorch - CPU or CUDA depending on TORCH_INDEX.
#   Pinned to REINVENT4's required version so the later `pip install -e .`
#   sees the requirement satisfied and keeps THIS build (CPU or CUDA) instead
#   of pulling its own.
RUN conda run -n reinvent4 pip install torch==2.12.0 torchvision \
    --index-url ${TORCH_INDEX}

# Layer 4: Clone REINVENT4
RUN git clone https://github.com/MolecularAI/REINVENT4.git /opt/reinvent4

# Layer 5: Install REINVENT4 and report which torch build is active
WORKDIR /opt/reinvent4
RUN conda run -n reinvent4 pip install -e . && \
    conda run -n reinvent4 python -c "import reinvent; print('REINVENT4 installed successfully')" && \
    conda run -n reinvent4 python -c "import torch; print('torch', torch.__version__, '| CUDA build:', torch.version.cuda or 'CPU-only')"

# Layer 6: AutoDock Vina 1.1.2 (bioconda 'autodock-vina' == 1.1.2) + chem helpers
RUN conda run -n reinvent4 mamba install -c bioconda autodock-vina -y && \
    conda run -n reinvent4 pip install rdkit scipy pandas

# Layer 7: DockStream (separate package + its own Python 3.7 env named "DockStream")
RUN git clone https://github.com/MolecularAI/DockStream.git /opt/DockStream
RUN mamba env create -f /opt/DockStream/environment.yml && \
    /opt/conda/envs/DockStream/bin/python /opt/DockStream/docker.py -h >/dev/null 2>&1 \
        && echo "DockStream entry point OK" || echo "WARN: docker.py -h returned nonzero (check env)"

# Layer 8: Workspace
WORKDIR /workspace
ENV PYTHONUNBUFFERED=1
RUN echo "source activate reinvent4" > ~/.bashrc
CMD ["/bin/bash"]

# Layer 9: Patch BucketCounter so the diversity filter survives checkpoint resume.
COPY scripts/patch_bucketcounter.py /tmp/patch_bucketcounter.py
RUN python /tmp/patch_bucketcounter.py
