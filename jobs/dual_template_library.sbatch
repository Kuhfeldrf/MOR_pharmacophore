#!/usr/bin/env bash
#SBATCH -J pharm_libB
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 120
#SBATCH -o logs/libB/%A_%a.log
#SBATCH --array=0-145%48
# Dual-template screen, stage 2: the full 14,595-target library under template B.
#
# Why this is a separate job. The E4' held-out percentile (design doc 6, T5) asks where
# DAMGO falls in the library score distribution, which only exists once the library has
# been scored against the same template. That is ~52 CPU-hours per seed - two orders of
# magnitude more than everything in jobs/dual_template.sbatch combined.
#
# Run jobs/dual_template.sbatch first. If T3 shows the two templates' dietary score
# vectors correlate at rho >= 0.9, template B is measuring what template A already
# measured, the dual-template design is unjustified, and this job should not be run at
# all. Spending the library screen before checking that is the expensive mistake.
#
# Mirrors jobs/screen_array.sbatch exactly except for template and output, so the two
# distributions are comparable. Window 8 and 800 conformers to match the dietary runs;
# note the original library screen used --n-conformers 50 at window 5, so DO NOT compare
# percentiles across those two settings (D4.2: ensemble size shifts the distribution).
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/jobs/env.sh"
ROOT="$MOR_ROOT"
cd "$ROOT"
mkdir -p logs/libB
python scripts/04_screen.py \
    --targets data/targets.tsv \
    --template "${TEMPLATE:-results/template_8f7r_p4/mor_template_excl_a1.phar}" \
    --with-exclusion \
    --block "${SLURM_ARRAY_TASK_ID}" \
    --nblocks "${NBLOCKS:-146}" \
    --out "${OUTDIR:-results/library_b8excl}" \
    --window "${WINDOW:-8}" \
    --n-conformers "${NCONF:-800}" \
    --seed "${SEED:-61453}"
