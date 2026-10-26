#!/bin/bash
#SBATCH --job-name=mor_resel
#SBATCH --partition=short
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=120G
#SBATCH --output=logs/report/%j.log
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/jobs/env.sh"
python scripts/16_reselect_with_e6.py --library results/library_w8 --bench-a results/stageB/nc800_w8_s20261020 --bench-b results/stageB/nc800_w8_s61453 --out docs
echo RESEL_DONE
