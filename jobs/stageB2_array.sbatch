#!/usr/bin/env bash
#SBATCH -J pharmB2
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH -t 240
#SBATCH -o logs/stageB2/%A_%a.log
#SBATCH --array=0-23%24
# Window 8 at 800 conformers, both seeds, full benchmark.
#
# Window 8 is the configuration that answers the study's central test, but its
# seed-to-seed Spearman on the full benchmark is 0.945 - just under the pre-registered
# 0.95 gate, which it cleared (0.960) on the shorter selection set. Rather than relax
# the gate, double the ensemble and re-measure.
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/jobs/env.sh"
ROOT="$MOR_ROOT"
cd "$ROOT"; mkdir -p logs/stageB2
SEEDS=(61453 20261020); NBLOCKS=12
i=${SLURM_ARRAY_TASK_ID}
blk=$(( i % NBLOCKS )); s_i=$(( i / NBLOCKS ))
SEED=${SEEDS[$s_i]}
python scripts/04_screen.py \
    --targets data/benchmark/benchmark.tsv \
    --template results/01_template/mor_template.phar \
    --block "$blk" --nblocks "$NBLOCKS" \
    --out "results/stageB/nc800_w8_s${SEED}" \
    --n-conformers 800 --window 8 --seed "$SEED"
