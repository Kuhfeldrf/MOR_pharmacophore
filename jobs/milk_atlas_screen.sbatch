#!/bin/bash
#SBATCH --job-name=mor_milkatlas
#SBATCH --partition=normal
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=logs/milk_atlas/%A_%a.log
#SBATCH --array=0-7%4
# Screen the milk opioid-motif atlas against both endogenous templates.
#
# Targets: data/milk_atlas/milk_atlas_targets.tsv - peptides that (a) carry a tier-1
# opioid motif, (b) were OBSERVED by MS in infant digesta, and (c) carry the motif at
# their own N-terminus, so the anchoring Tyr has a free alpha-amine. Built by
# scripts/25_milk_motif_atlas.py; rationale in docs/milk_atlas.md.
#
# This set contains its own positive controls: beta-casomorphin-4/5, neocasomorphin-6,
# alpha-lactorphin and the alphas1-casein exorphin all appear as released N-termini and
# all have published activity. A screen that fails to rank them above the novel sites is
# telling you about the screen, not about the novel sites. Check them first.
#
# Requires templates built by jobs/dual_template.sbatch (template B) and the existing
# results/template_p8 (template A).
#
# --- Why this is an array, in `normal`, and not one job in `short` -----------------
# It was one sequential job of 8 screens in `short`, and job 242432 was killed at the 4 h
# wall clock having finished NONE of them: 1.5 h to build the 35-peptide conformer
# ensemble, then 2.5 h of scoring that got 28 of 35 peptides into scores.tab before the
# clock ran out. One arm is therefore ~4.6 h and the whole loop ~32 h, so no arrangement
# of this script could ever have fitted in `short`, and the single-job form threw away
# 4 h of compute for zero completed output. Same wall-clock trap as the one described at
# the top of jobs/library_screen.sbatch, which is the second time on this project.
#
# One arm per array task: 8 tasks = 4 arms x 2 seeds, each with the full 10 h to itself
# (2.2x the measured arm, because the bare arms are unmeasured and the aromatic-rich
# candidates are the slow ones). `normal` is 24 h; `short` is 4 h and is not an option.
#
# Idempotent at ARM granularity, so a partial campaign resumes instead of restarting. An
# arm is finished only when scripts/04_screen.py has written block0000.json, which it does
# last (after the TSV, after scoring); a bare block0000/ directory means an interrupted
# arm and is redone from the start.
#
# There is no resume WITHIN an arm, and the conformer ensemble is not cached: the
# existence check at 04_screen.py:86 is obabel's "produced nothing" error check inside
# protonate_sdf, not a reuse cache, and build_block_sdf rebuilds conformers.sdf
# unconditionally. Measured here as 5431s of the interrupted arm's wall clock, spent
# rebuilding 666 MB that was already on disk.
#
# That is also the obvious saving, and it is deliberately not taken here: conformers
# depend on (targets, seed, window, n_conformers) and NOT on the template, so all four
# arms of a seed rebuild the same ensemble and the campaign pays for it 4x per seed,
# ~12 h across the eight arms. A cache would cost 666 MB per seed kept on a scratch
# filesystem that is 96% full and whose SDFs are transient by design (--keep-sdf,
# "scratch is tight"). Changing that is a change to 04_screen.py and affects every
# screen in the project, so it belongs in its own commit with its own measurement.
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/jobs/env.sh"
mkdir -p logs/milk_atlas

TA=results/template_p8
TB=results/template_8f7r_p4
SEEDS=(20261020 61453)
TARGETS=data/milk_atlas/milk_atlas_targets.tsv

TASK="${SLURM_ARRAY_TASK_ID:?this script must be submitted as an array}"
SEED="${SEEDS[$((TASK % 2))]}"

# Arm 0/1 carry exclusion volumes; arm 2/3 are the SAME targets WITHOUT them.
#
# The bare arms are not redundant and not optional. Two of the five novel candidates are
# aromatic-rich (YPYAVWYYPQ and MUC1 YPTYHTHGRYVPP, 4 aromatics each) and so sit exactly
# where this project's documented confound lives: score correlates with aromatic count at
# rho = +0.909, and a His-rich repeat once outscored every endogenous ligand
# (docs/optimisation_study.md 7.2). Exclusion volumes are the fix (D2.4), and the
# evidence that they worked on the dietary set was with-against-without: E8 AUROC 0.917
# with, 0.444 without. Without the no-exclusion arm there is nothing to difference
# against, and a high score for an aromatic-rich peptide is uninterpretable.
case $((TASK / 2)) in
  0) TAG=a8excl; TEMPLATE="$TA/mor_template_excl_a1.phar"; EX="--with-exclusion" ;;  # 8F7Q, beta-endorphin, extent 8
  1) TAG=b8excl; TEMPLATE="$TB/mor_template_excl_a1.phar"; EX="--with-exclusion" ;;  # 8F7R, endomorphin-1, complete ligand
  2) TAG=a8bare; TEMPLATE="$TA/mor_template.phar";         EX="" ;;
  3) TAG=b8bare; TEMPLATE="$TB/mor_template.phar";         EX="" ;;
  *) echo "FATAL: array index $TASK is outside 0-7." >&2; exit 1 ;;
esac

OUT="results/milk_atlas/${TAG}_s${SEED}"
if [ -f "$OUT/block0000.json" ]; then
    echo "SKIP tag=$TAG seed=$SEED - $OUT/block0000.json already written"
    exit 0
fi

echo "ARM tag=$TAG seed=$SEED template=$TEMPLATE out=$OUT"
# -u because 242432's log held nothing but the kill notice: the progress lines were still
# in a buffer when SLURM killed it, so there was no record of how far it had got.
python -u scripts/04_screen.py \
    --targets "$TARGETS" \
    --template "$TEMPLATE" \
    --out "$OUT" \
    --window 8 --n-conformers 800 --seed "$SEED" $EX

echo "MILK_ATLAS_ARM_DONE tag=$TAG seed=$SEED"
