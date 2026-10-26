# MOR_pharmacophore

Pharmacophore motif matching for peptide-receptor screening. Given a receptor with a known
endogenous peptide ligand, find dietary peptides that reproduce the geometry the ligand uses
to bind, rather than matching amino acid sequence.

Part of the Peptide-Receptor Discovery Pipeline (see `MOR_peptide_binding` and
`oprm1-opes-activation`). Runs on CPU only; it places no demand on the GPU partition.

## Method

1. **Interactions (PLIP).** Read the hydrogen bonds, hydrophobic contacts, salt bridges and
   aromatic stacking between the receptor and its own ligand in an experimental structure.
2. **Feature query (RDKit).** Encode the peptide as pharmacophore features
   (BaseFeatures.fdef) with 3D positions.
3. **Superposition and scoring (Align-it).** Score how closely a candidate peptide's chemical
   groups reproduce the template positions.
4. **Positional rule.** A match must satisfy the ligand's positional requirement; for the
   opioid receptors, a free N-terminus at the anchoring Tyr. A digestion peptide supports a
   Tyr site only if its N-terminus sits at that Tyr.

## Status

The screen has been run end to end against the mu-opioid receptor on ORCA. Plan and
pre-registered acceptance criteria: [`PLAN.md`](PLAN.md). Result:
[`docs/screen_validation.md`](docs/screen_validation.md).

| Step | In this repo | Tested |
|---|---|---|
| PLIP interactions, mu-opioid receptor / met-enkephalin (8F7Q chain P) | yes | yes |
| RDKit feature encoding | yes | yes |
| Free N-terminus rule | yes | yes |
| Align-it build from source (v1.0.4) | `env/build_alignit.sbatch` | yes |
| Align-it scoring wrapper | yes, score-table layout verified against output | yes |
| Template extraction from the experimental complex | `scripts/03_template.py` | yes |
| Library ingestion (14,571 peptides) | `scripts/01_library.py` | yes |
| Control sets, frozen before scoring | `scripts/02_controls.py` | yes |
| Conformers + scoring, SLURM array | `scripts/04_screen.py`, `jobs/` | yes, full library |
| Validation against pre-registered criteria | `scripts/05_validate.py` | yes |
| Tiered scoring | not yet - decile tiers are not stable enough to use (C5) | - |
| Table S2 affinities, all four opioid receptors | `data/reference/wang2023_tableS2.csv` | yes, cross-checked |
| Dual-template screen (8F7Q + 8F7R), pre-registered | code, data and design: `docs/dual_template_design.md` | not yet scored |
| SMILES targets, for ligands with no one-letter code | `conformers.build_from_smiles` | yes |
| 83-receptor batch run, ARM build | not yet | planned for the seed-fund Q1-Q2 milestones |

### Which endogenous ligand should the template come from?

Five endogenous peptides bind the mu-opioid receptor, but only two have an experimentally
determined mOR-bound pose, so only two can be templates:

| Endogenous muOR ligand | muOR pKi | mOR-bound structure |
|---|---|---|
| Endomorphin-1 | **8.73** | **8F7R** |
| Leu-enkephalin | 7.97 | none |
| Met-enkephalin | 7.40 | none |
| beta-Endorphin (1-31) | 7.11 | **8F7Q** |
| Dynorphin A (1-13) | 6.02 | 8F7W (kappa, not mu) |

The pipeline's original template came from 8F7Q, which is the 4th-strongest of the five at
muOR and non-selective mu-vs-delta. A second template from 8F7R is built and
pre-registered but not yet scored; the argument, the predictions and what would falsify
them are in [`docs/dual_template_design.md`](docs/dual_template_design.md).

**First run, 2026-10-20: C1, C2 and C4 pass; C3 and C5 fail.** The template recovers an
independent experimental ligand (endomorphin-1, held-out structure 8F7R) at the ~99th
percentile in both seeds, and the positional rule separates peptides that start with the
opioid message from the 74 library peptides that merely span it. But the screen is **not**
shown to add discrimination beyond the free "starts with Tyr" baseline, and at 50
conformers the score is not reproducible enough to tier. Numbers, both seeds, and what
to fix first are in [`docs/screen_validation.md`](docs/screen_validation.md).

The sequence-level pre-screen of milk proteins (39 opioid-like motif sites, 11 released as
N-termini) lives in `milk-opioid-motif-atlas`.

## Run the tests

```
pip install -r requirements.txt
pytest
```

## Run the screen

CPU only; it places no demand on the GPU partition.

```bash
sbatch env/build_env.sbatch          # conda env (a login-node solve is OOM-killed)
sbatch env/build_alignit.sbatch      # Align-it from source; needs the conda toolchain
python scripts/01_library.py --xlsx <peptidomics.xlsx>
python scripts/02_controls.py
python scripts/03_template.py        # PLIP + Align-it template from 8F7Q
python scripts/03b_targets.py
sbatch jobs/screen_array.sbatch      # 146 blocks, ~52 CPU-hours for 14,595 peptides
python scripts/05_validate.py --screen2 <second-seed run>
```

Test data: `data/8f7q_complex.pdb` is receptor chain R plus peptide chain P (residues 1-5,
Tyr-Gly-Gly-Phe-Met) of PDB 8F7Q, copied from `met-enkephalin-pose-overlay`.

## Licenses

Align-it is GPL-3.0 and is run as an external program, not bundled. PLIP, RDKit and OpenBabel
have their own licenses.
