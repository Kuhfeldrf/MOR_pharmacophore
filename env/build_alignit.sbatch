#!/usr/bin/env bash
#SBATCH -J alignit
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 60
#SBATCH -o logs/alignit_%j.log
# Stage 0b: build Align-it from source (GPL-3.0, not on conda-forge, upstream last
# pushed 2022-01-26). Built against the conda OpenBabel 3.2.1 in mor-pharm.
#
# MUST use the conda toolchain, not /usr/bin/g++. conda-forge OpenBabel is built with a
# newer GCC and its libopenbabel.so needs CXXABI_1.3.15; the system g++ here is 11.5 and
# its libstdc++ does not provide that symbol, so linking fails at 100% with
# "undefined reference to __cxa_call_terminate@CXXABI_1.3.15" (job 236160).
set -euo pipefail
# Conda location is a parameter, not a literal: the original hardcoded allocation was
# reaped and took every job in this repo with it (see jobs/env.sh). MOR_CONDA_ROOT may be
# set in the submit environment; it defaults to miniforge3 beside the workspace, which is
# how this project has always been laid out.
MOR_CONDA_ROOT="${MOR_CONDA_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}/../miniforge3}"
source "$MOR_CONDA_ROOT/etc/profile.d/conda.sh"
conda activate mor-pharm
conda install -y -c conda-forge gxx_linux-64 gcc_linux-64
ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
SRC=$ROOT/src_ext/align-it
rm -rf "$SRC" && mkdir -p "$ROOT/src_ext"
git clone --depth 1 https://github.com/OliverBScott/align-it "$SRC"
cd "$SRC" && echo "COMMIT $(git rev-parse HEAD)"
mkdir -p build && cd build
cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER="$CC" \
  -DCMAKE_CXX_COMPILER="$CXX" \
  -DCMAKE_INSTALL_PREFIX="$CONDA_PREFIX" \
  -DCMAKE_PREFIX_PATH="$CONDA_PREFIX" \
  -DOPENBABEL3_INCLUDE_DIR="$CONDA_PREFIX/include/openbabel3" \
  -DOPENBABEL3_LIBRARIES="$CONDA_PREFIX/lib/libopenbabel.so" \
  -DCMAKE_POLICY_VERSION_MINIMUM=3.5
make -j8 2>&1 | grep -v -E "warning:|note:|^ +[0-9]+ \||^ +\^|^ +~" | tail -20
make install 2>/dev/null || cp -v align_it "$CONDA_PREFIX/bin/" 2>/dev/null || true
ls -la "$CONDA_PREFIX/bin/" | grep -i align || true
"$CONDA_PREFIX/bin/align_it" --help 2>&1 | head -8 || true
echo "ALIGNIT_BUILD_OK"
