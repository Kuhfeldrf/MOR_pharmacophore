#!/bin/bash
#SBATCH --job-name=mor_ensz
#SBATCH --partition=short
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --output=logs/report/%j.log
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/jobs/env.sh"
python scripts/21_ensemble_size_statistic.py --stage results/stageB --out docs
echo ENSZ_DONE
