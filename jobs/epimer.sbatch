#!/bin/bash
#SBATCH --job-name=mor_epimer
#SBATCH --partition=short
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=logs/epimer/%j.log
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/jobs/env.sh"
python scripts/04_screen.py --targets data/epimers/epimer_pairs.tsv \
  --template results/01_template/mor_template.phar \
  --out results/epimers --window 8 --n-conformers 800 --seed 61453
echo EPIMER_DONE
