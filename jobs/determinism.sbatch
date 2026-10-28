#!/bin/bash
#SBATCH --job-name=mor_determinism
#SBATCH --partition=short
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=logs/determinism/%j.log
# C1/C2 of scripts/31_determinism_check.py: is the screen bit-for-bit reproducible
# within one environment, given a fixed seed and a fixed input?
#
# Why this exists as a committed script. Job 242912 ran this measurement from an inline
# `sbatch --wrap`, which left `scontrol show job` reporting Command=(null) and nothing in
# the repository to re-run. A determinism result whose own invocation is unrecorded is
# the same failure it was written to investigate: the regeneration discrepancy in
# docs/dual_template_result.md could not be attributed because the run did not record
# what produced it.
#
# Why `short` is the right partition here, unlike for the screens. Three peptides of 4-8
# residues, two 800-conformer builds and two align-it passes each. Window 8 of YPSYGLNY
# is the longest molecule involved, against the 44-rotatable-bond peptides the dietary
# screen embeds, so this is minutes of work, not hours. The 4 h ceiling that made
# jobs/milk_atlas_screen.sbatch unrunnable is not a constraint at this size.
#
# Scoring uses the extent-8 template WITH exclusion volumes, matching the runs whose
# discrepancy prompted this (D2.4: a bare template is known-broken, AUROC 0.444 against
# 0.917, so a determinism result on one would not describe any screen this project runs).
#
# Idempotent on the output file, so a resubmission does not redo a completed measurement.
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/jobs/env.sh"
mkdir -p logs/determinism

OUT=docs/determinism_check.json
TEMPLATE=results/template_p8/mor_template_excl_a1.phar

if [ -f "$OUT" ]; then
    echo "SKIP - $OUT already written; delete it to re-measure"
    exit 0
fi
if [ ! -f "$TEMPLATE" ]; then
    echo "FATAL: $TEMPLATE not found. Build it with jobs/dual_template.sbatch step 0." >&2
    exit 1
fi

# `python -u` because job 242432's log held only its kill notice: a buffered stdout is
# lost when a job is killed, and the per-peptide C1/C2 lines are the only running account
# of a job that can otherwise look stalled.
python -u scripts/31_determinism_check.py --template "$TEMPLATE" --out "$OUT"

echo DETERMINISM_DONE
