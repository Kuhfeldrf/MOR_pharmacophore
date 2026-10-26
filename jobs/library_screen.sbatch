#!/usr/bin/env bash
#SBATCH -J pharm_lib
#SBATCH -p normal
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH -t 1200
#SBATCH -o logs/library/%A_%a.log
#SBATCH --array=0-299%48
# Full library screen at the VALIDATED configuration (docs/optimisation_study.md):
#
#   window 8, 800 conformers, TVERSKY_REF, mean-of-ensemble
#
# Chosen because it is the only configuration that both passes the pre-registered
# reproducibility gate (seed-to-seed Spearman 0.983) and can be evaluated on the
# potency-pair control, which it passes (AUROC 0.803, CI [0.561, 0.970]).
#
# ~4,700 core-hours over 14,595 targets, 300 blocks of ~49 peptides.
#
# Runs in `normal` (24 h), not `short` (4 h). Measured on the benchmark this
# configuration costs ~145 s per peptide on 8 cores, so a 49-peptide block is ~2 h - but
# library peptides have median length 15 against the benchmark's ~7, so nearly all of
# them fill the whole 8-residue window and cost more than the benchmark average. The
# estimate lands close enough to 4 h that `short` would risk losing whole blocks to the
# wall clock, which is how this project lost a campaign once already.
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/jobs/env.sh"
ROOT="$MOR_ROOT"
cd "$ROOT"; mkdir -p logs/library
python scripts/04_screen.py \
    --targets data/targets.tsv \
    --template results/01_template/mor_template.phar \
    --block "${SLURM_ARRAY_TASK_ID}" --nblocks 300 \
    --out "${OUTDIR:-results/library_w8}" \
    --n-conformers 800 --window 8 --seed "${SEED:-61453}"
