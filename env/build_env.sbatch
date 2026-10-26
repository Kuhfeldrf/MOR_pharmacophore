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
    python=3.11 rdkit=2024.09.6 openbabel=3.2.1 cmake make pytest numpy pandas openpyxl \
    lxml
conda activate mor-pharm
# PLIP is noarch on conda-forge but pip keeps it off the solver critical path.
#
# --no-deps is load-bearing. A plain `pip install plip` resolves its openbabel
# requirement to the PyPI `openbabel` wheel, which bundles its own libopenbabel and its
# own 128 format plugins built with a different C++ ABI, and installs them *over* the
# conda-forge Python bindings in this same prefix. align_it is linked against the
# conda library, the wheel sets BABEL_LIBDIR to its own plugin directory on import, and
# every align_it call that reads a structure then dies on an undefined symbol -- which
# is how job 242423 failed after this env reported ENV_BUILD_OK.
#
# PLIP's actual needs are openbabel (conda-forge, above), lxml (conda-forge, above) and
# numpy (above), so --no-deps leaves nothing missing. src/mor_pharmacophore/alignit.py
# also pins BABEL_LIBDIR defensively; both guards are deliberate.
pip install --no-deps "plip==3.0.1"
# Fail here rather than three jobs later if the wheel ever sneaks back in.
if python -c "import openbabel,sys; sys.exit(0 if 'site-packages/openbabel' in openbabel.__file__ and __import__('pathlib').Path(openbabel.__file__).parent.joinpath('lib','libopenbabel.so').exists() else 1)" 2>/dev/null; then
    echo "FATAL: the pip openbabel wheel is installed and will break align_it." >&2
    echo "  Remove it and reinstall conda-forge openbabel:" >&2
    echo "    pip uninstall -y openbabel && conda install -y --force-reinstall -c conda-forge openbabel=3.2.1" >&2
    exit 1
fi
python -c "import rdkit, plip, openbabel; print(\"rdkit\", rdkit.__version__); print(\"plip ok\")"
obabel -V
cmake --version | head -1
echo "ENV_BUILD_OK"
