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

| Step | In this repo | Tested |
|---|---|---|
| PLIP interactions, mu-opioid receptor / met-enkephalin (8F7Q chain P) | yes | yes |
| RDKit feature encoding | yes | yes |
| Free N-terminus rule | yes | yes |
| Align-it scoring | wrapper stub only | skipped unless `align-it` is installed |
| Tiered scoring, 83-receptor batch run, ARM build | not yet | planned for the seed-fund Q1-Q2 milestones |

The sequence-level pre-screen of milk proteins (39 opioid-like motif sites, 11 released as
N-termini) lives in `milk-opioid-motif-atlas`.

## Run the tests

```
pip install -r requirements.txt
pytest
```

Test data: `data/8f7q_complex.pdb` is receptor chain R plus peptide chain P (residues 1-5,
Tyr-Gly-Gly-Phe-Met) of PDB 8F7Q, copied from `met-enkephalin-pose-overlay`.

## Licenses

Align-it is GPL-3.0 and is run as an external program, not bundled. PLIP, RDKit and OpenBabel
have their own licenses.
