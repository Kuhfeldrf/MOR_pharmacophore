#!/bin/bash
#SBATCH --job-name=mor_libreport
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --output=logs/report/%j.log
#
# Two steps, in order:
#   1. Re-run the evaluation over every stage-B configuration, so sweep_results.json
#      holds all four (including window 8 / 800 conformers) rather than the three a
#      later partial run left behind. Selection is still read from
#      docs/selection_results.json, so only the pre-selected metric/statistic per
#      window is scored.
#   2. Rank and tier the full library. Reads ~1.2 GB of per-conformer score tables,
#      which is why this is a batch job and not a login-node command.
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/jobs/env.sh"

python scripts/12_analyze_sweep.py --sweep results/stageB --out docs

python scripts/13_library_report.py \
  --library results/library_w8 \
  --benchmark-run results/stageB/nc800_w8_s61453 \
  --out docs

echo REPORT_DONE
