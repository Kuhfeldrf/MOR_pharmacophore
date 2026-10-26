#!/usr/bin/env bash
#SBATCH -J pharm_env
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 180
#SBATCH -o logs/build_env_%j.log
# Stage 0a: the mor-pharm environment. A login-node solve is OOM-killed on this
# cluster (see MOR_peptide_binding README), so this is a job by necessity.
set -euo pipefail
# Conda location is a parameter, not a literal: the original hardcoded allocation was
# reaped and took every job in this repo with it (see jobs/env.sh). MOR_CONDA_ROOT may be
# set in the submit environment; it defaults to miniforge3 beside the workspace, which is
# how this project has always been laid out.
MOR_CONDA_ROOT="${MOR_CONDA_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}/../miniforge3}"
source "$MOR_CONDA_ROOT/etc/profile.d/conda.sh"
conda create -y -n mor-pharm -c conda-forge \
    python=3.11 rdkit=2024.09.6 openbabel=3.2.1 cmake make pytest numpy pandas openpyxl
conda activate mor-pharm
# PLIP is noarch on conda-forge but pip keeps it off the solver critical path.
pip install "plip==3.0.1"
python -c "import rdkit, plip, openbabel; print(\"rdkit\", rdkit.__version__); print(\"plip ok\")"
obabel -V
cmake --version | head -1
echo "ENV_BUILD_OK"
