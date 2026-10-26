#!/bin/bash
#SBATCH --job-name=mor_stageC
#SBATCH --partition=normal
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --array=0-31
#SBATCH --output=logs/stageC/%A_%a.log
# Stage C: the full benchmark under the CORRECTED template, with and without
# receptor exclusion volumes, two seeds each, 8 blocks.
#
# The 14-peptide diagnostic showed exclusion volumes invert the message-vs-His-rich
# ordering (gap -0.0444 -> +0.0677). This run is what decides whether that survives
# on the pre-registered benchmark: S1, E1, E3 and E6 all need the full set, and the
# seed jitter has to be RE-MEASURED under exclusion volumes, because the scores
# compress and the old 0.0013 no longer applies.
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/jobs/env.sh"
NBLOCKS=8
BLOCK=$((SLURM_ARRAY_TASK_ID % NBLOCKS))
CFG=$((SLURM_ARRAY_TASK_ID / NBLOCKS))
SEEDS=(20261020 61453)
case $CFG in
  0) T=results/template_p8/mor_template.phar;         TAG=p8;     EX="";                SEED=${SEEDS[0]} ;;
  1) T=results/template_p8/mor_template.phar;         TAG=p8;     EX="";                SEED=${SEEDS[1]} ;;
  2) T=results/template_p8/mor_template_excl_a1.phar; TAG=p8excl; EX="--with-exclusion"; SEED=${SEEDS[0]} ;;
  3) T=results/template_p8/mor_template_excl_a1.phar; TAG=p8excl; EX="--with-exclusion"; SEED=${SEEDS[1]} ;;
esac
python scripts/04_screen.py --targets data/benchmark/benchmark.tsv --template $T \
  --out results/stageC/${TAG}_w8_s${SEED} \
  --window 8 --n-conformers 800 --seed $SEED --block $BLOCK --nblocks $NBLOCKS $EX
echo STAGEC_DONE cfg=$TAG seed=$SEED block=$BLOCK
