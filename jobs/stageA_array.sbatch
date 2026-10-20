#!/usr/bin/env bash
#SBATCH -J pharmA
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH -t 240
#SBATCH -o /scratch/kuhfeldr-Kuhfeld_temp/MOR_pharmacophore/logs/stageA/%A_%a.log
#SBATCH --array=0-71%24
# STUDY.md 4, stage A: the selection sweep.
#
# Both selection metrics (S1 reproducibility, S2 actives-vs-scrambles) are computable
# from the 72-peptide selection set alone, so the expensive axes are swept here and the
# full benchmark is run only at the configuration this chooses. The evaluation metrics
# are deliberately NOT computable from this set.
#
# 3 windows x 3 conformer counts x 2 seeds x 4 blocks = 72 tasks.
set -euo pipefail
source /scratch/kuhfeldr-Kuhfeld_temp/miniforge3/etc/profile.d/conda.sh
conda activate mor-pharm
ROOT=/scratch/kuhfeldr-Kuhfeld_temp/MOR_pharmacophore
cd "$ROOT"
mkdir -p logs/stageA
# Explicit: RDKit numThreads=0 would mean all 64 hardware cores and oversubscribe the
# node (measured CPULoad 319 in sweep 237134).
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}

WINDOWS=(5 8 12)
NCONFS=(50 150 400)
SEEDS=(61453 20261020)
NBLOCKS=4

i=${SLURM_ARRAY_TASK_ID}
blk=$(( i % NBLOCKS )); r1=$(( i / NBLOCKS ))
s_i=$(( r1 % 2 ));      r2=$(( r1 / 2 ))
nc_i=$(( r2 % 3 ));     w_i=$(( r2 / 3 ))

W=${WINDOWS[$w_i]}; NC=${NCONFS[$nc_i]}; SEED=${SEEDS[$s_i]}
TAG="nc${NC}_w${W}_s${SEED}"
echo "stageA task $i: window=$W n_conformers=$NC seed=$SEED block=$blk -> $TAG"

python scripts/04_screen.py \
    --targets data/benchmark/selection_set.tsv \
    --template results/01_template/mor_template.phar \
    --block "$blk" --nblocks "$NBLOCKS" \
    --out "results/stageA/$TAG" \
    --n-conformers "$NC" --window "$W" --seed "$SEED"
