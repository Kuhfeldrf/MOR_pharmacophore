#!/usr/bin/env bash
# Submit the whole dual-template campaign as one dependent chain. Run from the repo root
# on ORCA (a login node is fine; this only submits).
#
#   bash jobs/run_all_screens.sh
#
# Order matters, and the ordering is the point: everything cheap runs first, and the one
# expensive job is submitted by a script that reads a number rather than by a human who
# has stopped reading.
#
#   1  dual_template          build template B + exclusion volumes; dietary B8/B4 and
#                             held-out ligands, both seeds.                  ~minutes
#   2  analyse                scripts/24 (consensus, T1-T4, T6) and scripts/23
#                             (single-template E8/E9 on the new runs), then the T3 gate.
#                             If the gate passes it submits step 4 itself.    ~seconds
#   3  milk_atlas             35 released motif-bearing peptides vs both templates.
#                             Independent of the gate - it is cheap and it carries its
#                             own positive controls.                         ~minutes
#   4  library (GATED)        all 14,595 experimentally observed targets under template
#                             B, for the E4' percentile only.           ~52 CPU-h/seed
#
# Step 4 is submitted ONLY if T3 shows the templates are non-redundant. See
# scripts/26_t3_gate.py for why, and docs/dual_template_design.md 6 for the prediction.
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p logs/dual logs/milk_atlas logs/libB

J1=$(sbatch --parsable jobs/dual_template.sbatch)
echo "1  dual_template   $J1"

# Analysis + gate. Submitted as a job so it inherits the dependency; it is seconds of CPU.
J2=$(sbatch --parsable --dependency=afterok:"$J1" \
     --job-name=mor_dualanalyse --partition=short --time=00:20:00 \
     --cpus-per-task=2 --mem=8G --output=logs/dual/analyse_%j.log \
     --wrap "source \"\$SLURM_SUBMIT_DIR/jobs/env.sh\" && \
             python scripts/23_dietary_screen.py --runs results/dietary --out docs && \
             python scripts/24_dual_template_consensus.py --runs results/dietary \
                 --template-a p8excl --template-b b8excl \
                 --holdout results/holdout --out docs && \
             if python scripts/26_t3_gate.py --result docs/dual_template_result.json; then \
                 sbatch jobs/dual_template_library.sbatch; \
             else \
                 echo 'LIBRARY SCREEN NOT SUBMITTED - T3 gate failed'; \
             fi")
echo "2  analyse+gate    $J2   (submits step 4 if the gate passes)"

J3=$(sbatch --parsable --dependency=afterok:"$J1" jobs/milk_atlas_screen.sbatch)
echo "3  milk_atlas      $J3"

J3b=$(sbatch --parsable --dependency=afterok:"$J3" \
      --job-name=mor_atlasanalyse --partition=short --time=00:20:00 \
      --cpus-per-task=2 --mem=8G --output=logs/milk_atlas/analyse_%j.log \
      --wrap "source \"\$SLURM_SUBMIT_DIR/jobs/env.sh\" && \
              python scripts/27_milk_atlas_analysis.py --runs results/milk_atlas \
                  --targets data/milk_atlas/milk_atlas_targets.tsv --out docs")
echo "3b atlas analysis  $J3b"

cat <<EOF

Submitted. Watch with:  squeue -u \$USER
Results land in:
  docs/dietary_screen.json          single-template E8/E9, now including template B
  docs/dual_template_result.json    T1-T4, T6, consensus, E4'
  docs/milk_atlas_result.json       M1-M4 on the 35 motif-bearing released peptides
  results/milk_atlas/               35 motif-bearing released peptides, both templates
  results/library_b8excl/           only if the T3 gate passed

Read docs/dual_template_design.md 6 BEFORE the numbers, not after. Every prediction is
written down there already.
EOF
