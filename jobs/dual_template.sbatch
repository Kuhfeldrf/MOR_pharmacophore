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
source /scratch/kuhfeldr-Kuhfeld_temp/miniforge3/etc/profile.d/conda.sh
conda activate mor-pharm
cd "$SLURM_SUBMIT_DIR"
mkdir -p logs/dual
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

TB=results/template_8f7r_p4
SEEDS=(20261020 61453)

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
# B8 is primary: it matches A8 (results/dietary/p8excl_s*) on window, conformer count,
# seed, metric and exclusion treatment, so the only difference is the template.
# B4 is the extent-matched sensitivity run and is pre-registered as UNINTERPRETABLE as
# evidence for template B - at window 4 the aromatic confound favours the hypothesis and
# 8 of 20 peptides collapse onto a shared window. See design doc 5.2.
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
