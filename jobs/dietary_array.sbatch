#!/bin/bash
#SBATCH --job-name=mor_diet
#SBATCH --partition=short
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --array=0-5
#SBATCH --output=logs/dietary/%A_%a.log
# Dietary opioid peptides, 20 canonical sequences with measured IC50.
# Three templates x two seeds: the corrected 8-residue template is the one the
# structure supports (ARG213/ECL2 salt bridge present), p5 is the old truncated
# template kept for comparison, and p8+excl adds receptor exclusion volumes.
set -euo pipefail
source /scratch/kuhfeldr-Kuhfeld_temp/miniforge3/etc/profile.d/conda.sh
conda activate mor-pharm
cd "$SLURM_SUBMIT_DIR"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
CFG=$((SLURM_ARRAY_TASK_ID / 2))
SEEDS=(20261020 61453)
SEED=${SEEDS[$((SLURM_ARRAY_TASK_ID % 2))]}
case $CFG in
  0) T=results/template_p5/mor_template.phar;            TAG=p5;      EX="" ;;
  1) T=results/template_p8/mor_template.phar;            TAG=p8;      EX="" ;;
  2) T=results/template_p8/mor_template_excl_a1.phar;    TAG=p8excl;  EX="--with-exclusion" ;;
esac
python scripts/04_screen.py --targets data/dietary_targets.tsv --template $T \
  --out results/dietary/${TAG}_s${SEED} --window 8 --n-conformers 800 --seed $SEED $EX
echo DIET_DONE tag=$TAG seed=$SEED
