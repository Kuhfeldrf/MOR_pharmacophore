#!/bin/bash
#SBATCH --job-name=mor_dualtmpl
#SBATCH --partition=short
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=logs/dual/%j.log
# Dual-template screen, stage 1: build template B and run the dietary + held-out screens.
#
# Design and pre-registered predictions: docs/dual_template_design.md. Read that first -
# the window choice here is not free, and the reason it is 8 and not 4 is measured there.
#
# Template B comes from 8F7R (mOR - endomorphin-1 - Gi). Unlike 8F7Q there is no peptide
# extent to choose: the deposition resolves all four residues plus the C-terminal amide
# cap, so --peptide-residues 0 is the complete ligand, not a truncation.
#
# The library distribution needed for the E4' percentile (T5) is a separate, much longer
# job: jobs/dual_template_library.sbatch. Run this one first and check T1/T3 before
# spending 50+ CPU-hours on that.
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/jobs/env.sh"
mkdir -p logs/dual

TA=results/template_p8
TB=results/template_8f7r_p4
SEEDS=(20261020 61453)

# --- 0. Template A, rebuilt if absent --------------------------------------------
# Template A used to be assumed present, because it was. When the scratch allocation
# was reaped it went with everything else, and a clean clone has data/ but not results/
# (results/ is gitignored by D7.4). The holdout and milk-atlas screens both score
# against template A, so without this the campaign fails three jobs deep for a reason
# that has nothing to do with the campaign. Idempotent: skipped when already built.
if [ ! -f "$TA/mor_template_excl_a1.phar" ]; then
  echo "template A absent - rebuilding from data/8f7q_complex_p8.pdb"
  python scripts/03_template.py \
      --complex data/8f7q_complex_p8.pdb \
      --chain P \
      --out "$TA"
  python scripts/20_template_exclusion.py \
      --complex data/8f7q_complex_p8.pdb \
      --ligand-chain P \
      --template "$TA/mor_template.phar" \
      --out "$TA" \
      --alphas 0.5 1.0
else
  echo "template A already built, reusing $TA"
fi

# --- 1. Template B from the complete endomorphin-1 ligand -------------------------
# data/8f7r_complex_p4.pdb is built by scripts/22_prepare_complex.py and committed with
# its provenance JSON; rebuild it only if the deposition changes.
python scripts/03_template.py \
    --complex data/8f7r_complex_p4.pdb \
    --chain P \
    --out "$TB"

# --- 2. Exclusion volumes, mandatory ---------------------------------------------
# dietary_screen.md measured 0.917 with exclusion volumes against 0.444 without, same
# template and same set. A template without them is known-broken, so B is never run bare.
python scripts/20_template_exclusion.py \
    --complex data/8f7r_complex_p4.pdb \
    --ligand-chain P \
    --template "$TB/mor_template.phar" \
    --out "$TB" \
    --alphas 0.5 1.0

# --- 3. Dietary screen, both windows x both seeds ---------------------------------
# B8 is primary: it matches A8 (results/dietary/p8excl_s*, step 3a) on window, conformer
# count, seed, metric and exclusion treatment, so the only difference is the template.
# B4 is the extent-matched sensitivity run and is pre-registered as UNINTERPRETABLE as
# evidence for template B - at window 4 the aromatic confound favours the hypothesis and
# 8 of 20 peptides collapse onto a shared window. See design doc 5.2.

# --- 3a. Template A's dietary runs, rebuilt if absent -----------------------------
# Same reasoning as step 0, and the same bug: these were assumed present because they
# were, having been produced by jobs/dietary_array.sbatch CFG=2 in an earlier allocation.
# When that allocation was reaped the template was rebuilt by step 0 but these were not,
# so this job reported DUAL_TEMPLATE_DONE and the analysis behind it died on
# `no runs for template tag 'p8excl'` - three jobs deep, for a reason unrelated to the
# campaign, after every expensive screen had already run (job 242431).
#
# A8 is not optional: it is the control arm B8 is compared against, so T3 - and therefore
# the gate on the ~52 CPU-hour library screen - cannot be computed without it. Keeping it
# here rather than in the submitter means the comparison arm is built by the same job
# that builds the thing it is compared to, and cannot drift from it.
#
# Must stay identical to dietary_array.sbatch CFG=2: window 8, 800 conformers, alpha-1
# exclusion volumes, both seeds. Idempotent on block0000.json, which 04_screen.py writes
# last, so a completed run is never recomputed. ~1.9 h per seed when it does run.
for SEED in "${SEEDS[@]}"; do
  if [ -f "results/dietary/p8excl_s${SEED}/block0000.json" ]; then
    echo "template A dietary already screened for seed $SEED, reusing"
  else
    echo "template A dietary absent for seed $SEED - screening"
    python scripts/04_screen.py \
        --targets data/dietary/dietary_targets.tsv \
        --template "$TA/mor_template_excl_a1.phar" \
        --out "results/dietary/p8excl_s${SEED}" \
        --window 8 --n-conformers 800 --seed "$SEED" --with-exclusion
  fi
done

for SEED in "${SEEDS[@]}"; do
  for W in 8 4; do
    python scripts/04_screen.py \
        --targets data/dietary/dietary_targets.tsv \
        --template "$TB/mor_template_excl_a1.phar" \
        --out "results/dietary/b${W}excl_s${SEED}" \
        --window "$W" --n-conformers 800 --seed "$SEED" --with-exclusion
  done
done

# --- 4. Held-out ligands against both templates -----------------------------------
# DAMGO is built from SMILES (N-MePhe and glycinol have no one-letter code); every other
# target in this table is built from sequence. Same embedding protocol either way.
TA=results/template_p8
for SEED in "${SEEDS[@]}"; do
  python scripts/04_screen.py \
      --targets data/reference/holdout_targets.tsv \
      --template "$TA/mor_template_excl_a1.phar" \
      --out "results/holdout/a8excl_s${SEED}" \
      --window 8 --n-conformers 800 --seed "$SEED" --with-exclusion
  python scripts/04_screen.py \
      --targets data/reference/holdout_targets.tsv \
      --template "$TB/mor_template_excl_a1.phar" \
      --out "results/holdout/b8excl_s${SEED}" \
      --window 8 --n-conformers 800 --seed "$SEED" --with-exclusion
done

echo DUAL_TEMPLATE_DONE
