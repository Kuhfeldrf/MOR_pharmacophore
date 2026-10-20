#!/usr/bin/env bash
#SBATCH -J pharmB
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH -t 240
#SBATCH -o /scratch/kuhfeldr-Kuhfeld_temp/MOR_pharmacophore/logs/stageB/%A_%a.log
#SBATCH --array=0-71%24
# STUDY.md 4-5, stage B: evaluation.
#
# The three configurations chosen by the SELECTION metrics (one per window, STUDY.md
# 4.1) are run over the full 235-peptide benchmark, which adds the Tyr-matched decoys
# and the positional control that stage A deliberately excluded. Metric and ensemble
# statistic are read post hoc, so only (window, n_conformers, seed) vary here.
#
#   window  5, 400 conformers   chosen primary
#   window  8, 400 conformers   minimum window at which the potency pair is testable
#   window 12, 150 conformers   most family-pair comparisons
set -euo pipefail
source /scratch/kuhfeldr-Kuhfeld_temp/miniforge3/etc/profile.d/conda.sh
conda activate mor-pharm
ROOT=/scratch/kuhfeldr-Kuhfeld_temp/MOR_pharmacophore
cd "$ROOT"
mkdir -p logs/stageB
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}

WINDOWS=(5 8 12)
NCONFS=(400 400 150)
SEEDS=(61453 20261020)
NBLOCKS=12

i=${SLURM_ARRAY_TASK_ID}
blk=$(( i % NBLOCKS )); r1=$(( i / NBLOCKS ))
s_i=$(( r1 % 2 ));      cfg=$(( r1 / 2 ))

W=${WINDOWS[$cfg]}; NC=${NCONFS[$cfg]}; SEED=${SEEDS[$s_i]}
TAG="nc${NC}_w${W}_s${SEED}"
echo "stageB task $i: window=$W n_conformers=$NC seed=$SEED block=$blk -> $TAG"

python scripts/04_screen.py \
    --targets data/benchmark/benchmark.tsv \
    --template results/01_template/mor_template.phar \
    --block "$blk" --nblocks "$NBLOCKS" \
    --out "results/stageB/$TAG" \
    --n-conformers "$NC" --window "$W" --seed "$SEED"
