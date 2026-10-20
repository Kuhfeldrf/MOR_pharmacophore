#!/usr/bin/env bash
#SBATCH -J pharm_screen
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 120
#SBATCH -o /scratch/kuhfeldr-Kuhfeld_temp/MOR_pharmacophore/logs/screen/%A_%a.log
#SBATCH --array=0-145%48
# Stage 2+4: conformers and Align-it scoring, one block per array task.
#
# CPU only - this screen places no demand on the GPU partition, which is what the
# seed-fund narrative claims for it.
#
# 146 blocks x ~100 peptides covers the 14,595 targets. %48 matches the concurrency the
# sibling MOR_peptide_binding workflow uses on this cluster; blocks are independent, so
# raising it costs nothing but queue courtesy.
set -euo pipefail
source /scratch/kuhfeldr-Kuhfeld_temp/miniforge3/etc/profile.d/conda.sh
conda activate mor-pharm
ROOT=/scratch/kuhfeldr-Kuhfeld_temp/MOR_pharmacophore
cd "$ROOT"
mkdir -p logs/screen
# RDKit threading is set from the allocation, not left to default.
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
python scripts/04_screen.py \
    --targets data/targets.tsv \
    --template results/01_template/mor_template.phar \
    --block "${SLURM_ARRAY_TASK_ID}" \
    --nblocks "${NBLOCKS:-146}" \
    --out "${OUTDIR:-results/04_screen}" \
    --n-conformers "${NCONF:-50}" \
    --seed "${SEED:-61453}"
