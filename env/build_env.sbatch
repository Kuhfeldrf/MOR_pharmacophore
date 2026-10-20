#!/usr/bin/env bash
#SBATCH -J pharm_env
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 180
#SBATCH -o /scratch/kuhfeldr-Kuhfeld_temp/MOR_pharmacophore/logs/build_env_%j.log
# Stage 0a: the mor-pharm environment. A login-node solve is OOM-killed on this
# cluster (see MOR_peptide_binding README), so this is a job by necessity.
set -euo pipefail
source /scratch/kuhfeldr-Kuhfeld_temp/miniforge3/etc/profile.d/conda.sh
conda create -y -n mor-pharm -c conda-forge \
    python=3.11 rdkit=2024.09.6 openbabel=3.2.1 cmake make pytest numpy pandas openpyxl
conda activate mor-pharm
# PLIP is noarch on conda-forge but pip keeps it off the solver critical path.
pip install "plip==3.0.1"
python -c "import rdkit, plip, openbabel; print(\"rdkit\", rdkit.__version__); print(\"plip ok\")"
obabel -V
cmake --version | head -1
echo "ENV_BUILD_OK"
