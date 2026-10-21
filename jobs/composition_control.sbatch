#!/bin/bash
#SBATCH --job-name=mor_compctl
#SBATCH --partition=short
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=logs/report/%j.log
set -euo pipefail
source /scratch/kuhfeldr-Kuhfeld_temp/miniforge3/etc/profile.d/conda.sh
conda activate mor-pharm
cd "$SLURM_SUBMIT_DIR"
python scripts/14_composition_control.py --library results/library_w8 --benchmark-run results/stageB/nc800_w8_s61453 --out docs
echo CTL_DONE
