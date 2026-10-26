#!/bin/bash
#SBATCH --job-name=mor_milkatlas
#SBATCH --partition=short
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=logs/milk_atlas/%j.log
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
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/jobs/env.sh"
mkdir -p logs/milk_atlas

TA=results/template_p8
TB=results/template_8f7r_p4
SEEDS=(20261020 61453)

for SEED in "${SEEDS[@]}"; do
  # Template A: 8F7Q, beta-endorphin, extent 8
  python scripts/04_screen.py \
      --targets data/milk_atlas/milk_atlas_targets.tsv \
      --template "$TA/mor_template_excl_a1.phar" \
      --out "results/milk_atlas/a8excl_s${SEED}" \
      --window 8 --n-conformers 800 --seed "$SEED" --with-exclusion
  # Template B: 8F7R, endomorphin-1, complete ligand
  python scripts/04_screen.py \
      --targets data/milk_atlas/milk_atlas_targets.tsv \
      --template "$TB/mor_template_excl_a1.phar" \
      --out "results/milk_atlas/b8excl_s${SEED}" \
      --window 8 --n-conformers 800 --seed "$SEED" --with-exclusion

  # The SAME targets WITHOUT exclusion volumes, both templates.
  #
  # Not redundant, and not optional. Two of the five novel candidates are aromatic-rich
  # (YPYAVWYYPQ and MUC1 YPTYHTHGRYVPP, 4 aromatics each) and so sit exactly where this
  # project's documented confound lives: score correlates with aromatic count at
  # rho = +0.909, and a His-rich repeat once outscored every endogenous ligand
  # (docs/optimisation_study.md 7.2). Exclusion volumes are the fix (D2.4), and the
  # evidence that they worked on the dietary set was with-against-without: E8 AUROC 0.917
  # with, 0.444 without. Without the no-exclusion arm there is nothing to difference
  # against, and a high score for an aromatic-rich peptide is uninterpretable.
  python scripts/04_screen.py \
      --targets data/milk_atlas/milk_atlas_targets.tsv \
      --template "$TA/mor_template.phar" \
      --out "results/milk_atlas/a8bare_s${SEED}" \
      --window 8 --n-conformers 800 --seed "$SEED"
  python scripts/04_screen.py \
      --targets data/milk_atlas/milk_atlas_targets.tsv \
      --template "$TB/mor_template.phar" \
      --out "results/milk_atlas/b8bare_s${SEED}" \
      --window 8 --n-conformers 800 --seed "$SEED"
done

echo MILK_ATLAS_SCREEN_DONE
