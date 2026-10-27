#!/bin/bash
#SBATCH --job-name=mor_extent9
#SBATCH --partition=normal
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=logs/extent9/%j.log
# W1-W4 of docs/post_t3_design.md: the extent-9 template, and the dietary screen at
# window 9 that decides whether this design survives.
#
# Why extent 9 and not another template. Measured from the 8F7Q coordinates
# (scripts/29_ligand_contacts.py, docs/8f7q_contacts.tsv): every ECL2 and TM6-tip contact
# beta-endorphin makes is at residues 4-9 - Thr6/Arg213, Ser7/Lys305(6.58),
# Glu8/Arg213+Gln214, Lys9/Asp218 - and residues 10-21 touch only ECL1 and the TM1 tip,
# six of them nothing at all. Extent 8 already spans the Glu8 salt bridge; extent 9 adds
# Lys9-Asp218, the last one, and nothing beyond it is worth adding. See post_t3_design 2,
# which supersedes structural_basis 5 on this point.
#
# Why this is cheap and sequential rather than an array. Two seeds at ~2.5 h each plus a
# few minutes of template build fits inside one 10 h slot with room to spare, and a single
# job means the template cannot be built twice concurrently by two array tasks racing on
# the same output directory. `normal` because 2 x 2.5 h would not fit `short`'s 4 h - the
# mistake jobs/milk_atlas_screen.sbatch made and that cost job 242432 four hours for
# nothing.
#
# Idempotent per seed on block0000.json, which 04_screen.py writes last.
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/jobs/env.sh"
mkdir -p logs/extent9

T9=results/template_p9
SEEDS=(20261020 61453)

# --- Template, extent 9, with exclusion volumes ------------------------------------
# Exclusion volumes are not optional: dietary_screen.md measured AUROC 0.917 with them
# against 0.444 without on the same template and set, so a bare template is known-broken
# (D2.4). Rebuilt only when absent, so a resubmission does not redo it.
if [ ! -f "$T9/mor_template_excl_a1.phar" ]; then
    echo "building extent-9 template from data/8f7q_complex_p9.pdb"
    python -u scripts/03_template.py \
        --complex data/8f7q_complex_p9.pdb \
        --chain P \
        --out "$T9"
    python -u scripts/20_template_exclusion.py \
        --complex data/8f7q_complex_p9.pdb \
        --ligand-chain P \
        --template "$T9/mor_template.phar" \
        --out "$T9" \
        --alphas 0.5 1.0
else
    echo "extent-9 template already built, reusing $T9"
fi

# W2 is read off these two files by eye and by scripts/30_extent9_gates.py: the extent-9
# template must carry a feature the extent-8 one does not, or extent 9 is extent 8 with
# extra atoms.
echo "--- extent 9 template features ---"; grep -c . "$T9/mor_template.phar" || true
echo "--- extent 8 template features ---"; grep -c . results/template_p8/mor_template.phar || true

# --- Dietary screen at window 9, both seeds ----------------------------------------
# Window 9 to match the template extent, so reference and candidate describe the same
# amount of molecule (the extent-mismatch argument in scripts/22_prepare_complex.py).
for SEED in "${SEEDS[@]}"; do
    OUT="results/dietary/p9excl_s${SEED}"
    if [ -f "$OUT/block0000.json" ]; then
        echo "SKIP seed $SEED - $OUT/block0000.json already written"
        continue
    fi
    echo "screening dietary at window 9, seed $SEED"
    python -u scripts/04_screen.py \
        --targets data/dietary/dietary_targets.tsv \
        --template "$T9/mor_template_excl_a1.phar" \
        --out "$OUT" \
        --window 9 --n-conformers 800 --seed "$SEED" --with-exclusion
done

echo EXTENT9_DONE
