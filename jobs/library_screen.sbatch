#!/usr/bin/env bash
#SBATCH -J pharm_lib
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH -t 240
#SBATCH -o /scratch/kuhfeldr-Kuhfeld_temp/MOR_pharmacophore/logs/library/%A_%a.log
#SBATCH --array=0-299%48
# Full library screen at the VALIDATED configuration (docs/optimisation_study.md):
#
#   window 8, 800 conformers, TVERSKY_REF, mean-of-ensemble
#
# Chosen because it is the only configuration that both passes the pre-registered
# reproducibility gate (seed-to-seed Spearman 0.983) and can be evaluated on the
# potency-pair control, which it passes (AUROC 0.803, CI [0.561, 0.970]).
#
# ~4,700 core-hours over 14,595 targets; 300 blocks of ~49 peptides keeps each task
# inside the 4 h partition limit with margin.
set -euo pipefail
source /scratch/kuhfeldr-Kuhfeld_temp/miniforge3/etc/profile.d/conda.sh
conda activate mor-pharm
ROOT=/scratch/kuhfeldr-Kuhfeld_temp/MOR_pharmacophore
cd "$ROOT"; mkdir -p logs/library
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
python scripts/04_screen.py \
    --targets data/targets.tsv \
    --template results/01_template/mor_template.phar \
    --block "${SLURM_ARRAY_TASK_ID}" --nblocks 300 \
    --out "${OUTDIR:-results/library_w8}" \
    --n-conformers 800 --window 8 --seed "${SEED:-61453}"
