#!/bin/bash
#SBATCH --job-name=mor_shortlist
#SBATCH --partition=short
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=logs/report/%j.log
set -euo pipefail
source /scratch/kuhfeldr-Kuhfeld_temp/miniforge3/etc/profile.d/conda.sh
conda activate mor-pharm
cd "$SLURM_SUBMIT_DIR"
python scripts/15_shortlist.py --library results/library_w8 --out docs
echo SL_DONE
