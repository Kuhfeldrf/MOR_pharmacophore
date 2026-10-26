# Shared environment resolution for every SLURM job in this repo. Sourced, never executed.
#
# Why this file exists. Every job script used to carry
# `/scratch/kuhfeldr-Kuhfeld_temp/miniforge3/etc/profile.d/conda.sh` as a literal, in 24
# places across 23 files, plus the workspace root in 8 more and several absolute
# `#SBATCH --output` paths. When that scratch allocation was reaped the whole campaign
# became unrunnable, and because each script begins `set -euo pipefail` every job would
# have died on its first line with no useful message. A new allocation has a new name, so
# hardcoding it again just reschedules the same failure.
#
# Resolution order, all overridable from the submit environment:
#   MOR_ROOT       workspace root        default: $SLURM_SUBMIT_DIR, else $PWD
#   MOR_ENV        conda env name        default: mor-pharm
#   MOR_CONDA_SH   path to conda.sh      default: searched, see below
#
# Usage, after the #SBATCH block:
#   set -euo pipefail
#   source "${SLURM_SUBMIT_DIR:-$PWD}/jobs/env.sh"
#
# Keep `#SBATCH --output` paths RELATIVE. SLURM parses those directives before any shell
# runs, so a variable in them is not expanded; relative paths resolve against the
# submission directory, which is what we want anyway.

MOR_ROOT="${MOR_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
MOR_ENV="${MOR_ENV:-mor-pharm}"

if [ ! -d "$MOR_ROOT" ]; then
    echo "FATAL: MOR_ROOT=$MOR_ROOT does not exist." >&2
    exit 1
fi
cd "$MOR_ROOT"

# Locate conda. The sibling of the workspace is checked first because that is how this
# project has always been laid out (<alloc>/miniforge3 next to <alloc>/MOR_pharmacophore).
CONDA_SH="${MOR_CONDA_SH:-}"
if [ -z "$CONDA_SH" ]; then
    for candidate in \
        "$MOR_ROOT/../miniforge3" \
        "$MOR_ROOT/../miniconda3" \
        "$HOME/miniforge3" \
        "$HOME/miniconda3" \
        "$HOME/mambaforge"
    do
        if [ -f "$candidate/etc/profile.d/conda.sh" ]; then
            CONDA_SH="$candidate/etc/profile.d/conda.sh"
            break
        fi
    done
fi
if [ -z "$CONDA_SH" ] || [ ! -f "$CONDA_SH" ]; then
    echo "FATAL: no conda installation found." >&2
    echo "  Looked beside the workspace ($MOR_ROOT/..) and in \$HOME." >&2
    echo "  Set MOR_CONDA_SH=/path/to/conda.sh, or bootstrap with" >&2
    echo "  jobs/orca_bootstrap.sh which installs miniforge and builds the env." >&2
    exit 1
fi

# shellcheck disable=SC1090
source "$CONDA_SH"
if ! conda activate "$MOR_ENV" 2>/dev/null; then
    echo "FATAL: conda env '$MOR_ENV' not found in $(dirname "$(dirname "$CONDA_SH")")." >&2
    echo "  Available: $(conda env list 2>/dev/null | awk '!/^#/&&NF{print $1}' | tr '\n' ' ')" >&2
    echo "  Build it with env/build_env.sbatch, or jobs/orca_bootstrap.sh." >&2
    exit 1
fi

# RDKit's numThreads=0 means every core on the machine, not every core allocated, so the
# allocation is made explicit here for every job (see conformers.allocated_threads).
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

# Align-it is run as an external program; fail early and loudly rather than mid-screen.
if ! command -v align_it >/dev/null 2>&1; then
    echo "WARNING: align_it is not on PATH in env '$MOR_ENV'." >&2
    echo "  Scoring jobs will fail. Build it with env/build_alignit.sbatch." >&2
fi
