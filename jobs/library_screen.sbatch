#!/usr/bin/env bash
#SBATCH -J pharm_lib
#SBATCH -p normal
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH -t 1200
#SBATCH -o logs/library/%A_%a.log
#SBATCH --array=0-299%48
# Full library screen at the validated configuration:
#
#   window 8, 800 conformers, TVERSKY_REF, mean-of-ensemble,
#   extent-8 template WITH exclusion volumes (alpha 1.0)
#
# Window, ensemble size, metric and statistic are from docs/optimisation_study.md: the
# only configuration that both passes the pre-registered reproducibility gate
# (seed-to-seed Spearman 0.983) and can be evaluated on the potency-pair control, which
# it passes (AUROC 0.803, CI [0.561, 0.970]).
#
# ~4,700 core-hours over 14,595 targets, 300 blocks of ~49 peptides.
#
# Runs in `normal` (24 h), not `short` (4 h). Measured on the benchmark this
# configuration costs ~145 s per peptide on 8 cores, so a 49-peptide block is ~2 h - but
# library peptides have median length 15 against the benchmark's ~7, so nearly all of
# them fill the whole 8-residue window and cost more than the benchmark average. The
# estimate lands close enough to 4 h that `short` would risk losing whole blocks to the
# wall clock, which is how this project lost a campaign once already.
#
# --- Why the exclusion volumes are not optional, and why this script changed ----------
#
# Until 2026-10-28 this script screened `results/01_template/mor_template.phar` with no
# --with-exclusion flag, i.e. the bare template. Three independent results now say that
# configuration is broken, not merely suboptimal:
#
#   * docs/dietary_screen.md - AUROC 0.444 bare against 0.917 with exclusion volumes, on
#     the same template and the same 20 peptides (D2.4).
#   * docs/optimisation_study.md 7.2 - over all 14,571 library windows the bare score
#     correlates with aromatic count at rho = +0.909 and is close to a staircase in it.
#     Every endogenous ligand sits at 2-3 aromatics; every tier-1 window sits at 4-5,
#     in strata of 19 and 5 windows containing no known active.
#   * docs/milk_atlas_result.md - pre-registered M3, measured on 34 peptides: rho vs
#     aromatic count +0.805 bare against -0.012 with exclusion volumes, and the bare
#     arms rank the highest-aromatic novel windows first (MUC1, 4 aromatics, rank 1)
#     while the exclusion arms rank the positive controls first.
#
# So the ranking the bare configuration produced - generated 2026-10-21 - is a ranking of
# aromatic content. It must not be reused or extended, and rebuilding it is consequence 5
# of docs/dietary_screen.md 5. Those outputs are kept as the record of what the bare score
# does at library scale, and renamed docs/library_bare_* so that the name itself carries
# the warning: the report scripts write the unprefixed names, so leaving them in place
# would have let a rebuild silently overwrite the only surviving copy. The rebuild also
# writes its raw scores to a new directory, results/library_w8_excl.
#
# The old --template path is also gone: results/01_template went with the reaped scratch
# allocation and is not in version control. results/template_p8 is the extent-8 template
# every subsequent experiment in this project used, so the rebuild is on the same
# template as the dietary, dual-template and atlas results rather than on a third one.
#
# Idempotent per block on the block JSON, which 04_screen.py writes last, so a
# resubmission after a wall-clock loss redoes only the blocks that did not finish.
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/jobs/env.sh"
ROOT="$MOR_ROOT"
cd "$ROOT"; mkdir -p logs/library

TEMPLATE="${TEMPLATE:-results/template_p8/mor_template_excl_a1.phar}"
OUTDIR="${OUTDIR:-results/library_w8_excl}"
BLOCK=$(printf 'block%04d' "${SLURM_ARRAY_TASK_ID:?this script must be submitted as an array}")

if [ ! -f "$TEMPLATE" ]; then
    echo "FATAL: $TEMPLATE not found. Build it with jobs/dual_template.sbatch step 0." >&2
    exit 1
fi
if [ -f "$OUTDIR/$BLOCK.json" ]; then
    echo "SKIP $BLOCK - $OUTDIR/$BLOCK.json already written"
    exit 0
fi

# `python -u` because job 242432's log held only its kill notice: a buffered stdout is
# lost when a job is killed, and at ~2 h per block the progress lines are the only way to
# tell a slow block from a stalled one.
python -u scripts/04_screen.py \
    --targets data/targets.tsv \
    --template "$TEMPLATE" \
    --block "${SLURM_ARRAY_TASK_ID}" --nblocks 300 \
    --out "$OUTDIR" \
    --n-conformers 800 --window 8 --seed "${SEED:-61453}" \
    --with-exclusion
