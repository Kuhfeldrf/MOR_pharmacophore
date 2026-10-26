#!/usr/bin/env bash
#SBATCH -J pharm_sweep
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH -t 240
#SBATCH -o logs/sweep/%A_%a.log
#SBATCH --array=0-17%9
# STUDY.md 4: the optimisation sweep.
#
# 9 conformer configurations (n_conformers x window) x 2 seeds = 18 tasks, each over the
# whole 235-peptide benchmark in one block. Scoring metric and ensemble statistic are
# read post hoc from the retained per-conformer score tables, so they cost nothing here.
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/jobs/env.sh"
ROOT="$MOR_ROOT"
cd "$ROOT"
mkdir -p logs/sweep

NCONFS=(50 200 500)
WINDOWS=(5 8 0)
SEEDS=(61453 20261020)

i=${SLURM_ARRAY_TASK_ID}
seed_i=$(( i / 9 ))
rest=$(( i % 9 ))
nc_i=$(( rest / 3 ))
w_i=$(( rest % 3 ))

NC=${NCONFS[$nc_i]}
W=${WINDOWS[$w_i]}
SEED=${SEEDS[$seed_i]}
TAG="nc${NC}_w${W}_s${SEED}"

echo "sweep task $i: n_conformers=$NC window=$W seed=$SEED -> $TAG"
python scripts/04_screen.py \
    --targets data/benchmark/benchmark.tsv \
    --template results/01_template/mor_template.phar \
    --block 0 --nblocks 1 \
    --out "results/sweep/$TAG" \
    --n-conformers "$NC" \
    --window "$W" \
    --seed "$SEED"
