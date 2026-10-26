#!/bin/bash
#SBATCH --job-name=mor_diagexcl
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
T=results/01_template
# same conformers (same seed) scored against three templates
python scripts/04_screen.py --targets data/diagnostic/diag_targets.tsv --template $T/mor_template.phar \
  --out results/diag/base --window 8 --n-conformers 800 --seed 61453
python scripts/04_screen.py --targets data/diagnostic/diag_targets.tsv --template $T/mor_template_excl_a0p5.phar \
  --out results/diag/excl_a0p5 --window 8 --n-conformers 800 --seed 61453 --with-exclusion
python scripts/04_screen.py --targets data/diagnostic/diag_targets.tsv --template $T/mor_template_excl_a1.phar \
  --out results/diag/excl_a1 --window 8 --n-conformers 800 --seed 61453 --with-exclusion
echo DIAG_DONE
