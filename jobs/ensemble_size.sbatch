#!/bin/bash
#SBATCH --job-name=mor_ensz
#SBATCH --partition=short
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --output=logs/report/%j.log
set -euo pipefail
source /scratch/kuhfeldr-Kuhfeld_temp/miniforge3/etc/profile.d/conda.sh
conda activate mor-pharm
cd "$SLURM_SUBMIT_DIR"
python scripts/21_ensemble_size_statistic.py --stage results/stageB --out docs
echo ENSZ_DONE
