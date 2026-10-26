#!/usr/bin/env bash
#SBATCH -J pharm_screen
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 120
#SBATCH -o logs/screen/%A_%a.log
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
source "${SLURM_SUBMIT_DIR:-$PWD}/jobs/env.sh"
ROOT="$MOR_ROOT"
cd "$ROOT"
mkdir -p logs/screen
# RDKit threading is set from the allocation, not left to default.
python scripts/04_screen.py \
    --targets data/targets.tsv \
    --template results/01_template/mor_template.phar \
    --block "${SLURM_ARRAY_TASK_ID}" \
    --nblocks "${NBLOCKS:-146}" \
    --out "${OUTDIR:-results/04_screen}" \
    --n-conformers "${NCONF:-50}" \
    --seed "${SEED:-61453}"
