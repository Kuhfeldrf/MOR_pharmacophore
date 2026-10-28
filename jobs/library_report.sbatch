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
#
# --- What this now ranks, and what it will not re-rank -------------------------------
#
# This script used to read results/library_w8, the bare-template screen. That ranking is
# a ranking of aromatic content (rho = +0.909 over all 14,571 windows,
# docs/optimisation_study.md 7.2; see the header of jobs/library_screen.sbatch for the
# three results that establish it), so it must not be regenerated or extended. The
# defaults below point at the exclusion-volume rebuild instead.
#
# Both inputs are checked before any work starts, because both went with the reaped
# scratch allocation and neither is in version control. The rebuild chain is, in order:
#
#   1. a benchmark screen at window 8 / 800 conformers / exclusion volumes, which this
#      script needs for the tier edges and the decoy background (--benchmark-run);
#   2. jobs/library_screen.sbatch, 300 blocks, ~4,700 core-hours;
#   3. this job.
#
# Step 1 is not optional and is not included here: --benchmark-run is read
# unconditionally by scripts/13_library_report.py, and a tier edge drawn from a bare
# benchmark run against exclusion-volume library scores would be comparing two different
# score scales.
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/jobs/env.sh"

LIBRARY="${LIBRARY:-results/library_w8_excl}"
BENCHMARK_RUN="${BENCHMARK_RUN:-results/stageB/nc800_w8_s61453_excl}"

for d in "$LIBRARY" "$BENCHMARK_RUN"; do
    if [ ! -d "$d" ]; then
        echo "FATAL: $d not found." >&2
        echo "  Run the rebuild chain in this script's header first; neither input" >&2
        echo "  survived the reaped scratch allocation." >&2
        exit 1
    fi
done

python scripts/12_analyze_sweep.py --sweep results/stageB --out docs

python scripts/13_library_report.py \
  --library "$LIBRARY" \
  --benchmark-run "$BENCHMARK_RUN" \
  --out docs

echo REPORT_DONE
