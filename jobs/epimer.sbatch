#!/bin/bash
#SBATCH --job-name=mor_epimer
#SBATCH --partition=short
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=logs/epimer/%j.log
set -euo pipefail
source /scratch/kuhfeldr-Kuhfeld_temp/miniforge3/etc/profile.d/conda.sh
conda activate mor-pharm
cd "$SLURM_SUBMIT_DIR"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
python scripts/04_screen.py --targets data/epimer_targets.tsv \
  --template results/01_template/mor_template.phar \
  --out results/epimers --window 8 --n-conformers 800 --seed 61453
echo EPIMER_DONE
