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
source /scratch/kuhfeldr-Kuhfeld_temp/miniforge3/etc/profile.d/conda.sh
conda activate mor-pharm
cd "$SLURM_SUBMIT_DIR"
mkdir -p logs/milk_atlas
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

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
done

echo MILK_ATLAS_SCREEN_DONE
