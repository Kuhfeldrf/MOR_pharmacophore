#!/bin/bash
#SBATCH --job-name=mor_lowconf
#SBATCH --partition=short
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --array=0-7
#SBATCH --output=logs/lowconf/%A_%a.log
# Low conformer counts, the regime where McNutt et al. 2023 locate their
# larger-ensemble penalty (<25 vs hundreds). The 400-vs-800 test could not probe it.
set -euo pipefail
source /scratch/kuhfeldr-Kuhfeld_temp/miniforge3/etc/profile.d/conda.sh
conda activate mor-pharm
cd "$SLURM_SUBMIT_DIR"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
NCS=(10 25 50 100)
SEEDS=(20261020 61453)
NC=${NCS[$((SLURM_ARRAY_TASK_ID / 2))]}
SEED=${SEEDS[$((SLURM_ARRAY_TASK_ID % 2))]}
python scripts/04_screen.py --targets data/benchmark/benchmark.tsv \
  --template results/01_template/mor_template.phar \
  --out results/lowconf/nc${NC}_w8_s${SEED} \
  --window 8 --n-conformers $NC --seed $SEED
echo LOWCONF_DONE nc=$NC seed=$SEED
