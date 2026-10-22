#!/bin/bash
#SBATCH --job-name=mor_diagp8
#SBATCH --partition=short
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=logs/report/%j.log
set -euo pipefail
source /scratch/kuhfeldr-Kuhfeld_temp/miniforge3/etc/profile.d/conda.sh
conda activate mor-pharm
cd "$SLURM_SUBMIT_DIR"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
python scripts/04_screen.py --targets data/diag_targets.tsv \
  --template results/template_p8/mor_template.phar \
  --out results/diag/p8_base --window 8 --n-conformers 800 --seed 61453
python scripts/04_screen.py --targets data/diag_targets.tsv \
  --template results/template_p8/mor_template_excl_a1.phar \
  --out results/diag/p8_excl --window 8 --n-conformers 800 --seed 61453 --with-exclusion
python scripts/04_screen.py --targets data/diag_targets.tsv \
  --template results/template_p5/mor_template_excl_a1.phar \
  --out results/diag/p5_excl_newpocket --window 8 --n-conformers 800 --seed 61453 --with-exclusion
echo DIAGP8_DONE
