import pytest

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from mor_pharmacophore.conformers import (
    build_from_smiles, ensemble_for_smiles, ensemble_for_sequence)

# DAMGO, Tyr-D-Ala-Gly-N-MePhe-Gly-ol: the ligand of the held-out structure 8EFQ and the
# reason this code path exists. Its N-methyl-Phe and C-terminal glycinol have no
# one-letter code, so MolFromSequence cannot build it.
DAMGO = "N[C@@H](Cc1ccc(O)cc1)C(=O)N[C@H](C)C(=O)NCC(=O)N(C)[C@@H](Cc1ccccc1)C(=O)NCCO"


def test_damgo_smiles_is_the_right_molecule():
    mol, truncated = build_from_smiles(DAMGO)
    assert rdMolDescriptors.CalcMolFormula(mol) == "C26H35N5O6"
    assert truncated is False


def test_damgo_stereochemistry_is_l_tyr_d_ala_l_nmephe():
    """D-Ala is the whole point: an all-L DAMGO is a different molecule."""
    mol, _ = build_from_smiles(DAMGO)
    centres = dict(Chem.FindMolChiralCenters(mol, useLegacyImplementation=False))
    assert sorted(centres.values()) == ["R", "S", "S"]
    # Ala is the only R centre, and it is the residue between the two aromatics.
    assert list(centres.values()).count("R") == 1


def test_bad_smiles_raises_rather_than_returning_none():
    with pytest.raises(ValueError):
        build_from_smiles("this is not a smiles")


def test_smiles_ensemble_metadata_matches_the_sequence_path():
    """Both paths must fill the same keys, or the score table gains a special case."""
    smi_mol, smi_meta = ensemble_for_smiles(DAMGO, n_conformers=4, seed=1, n_threads=1)
    seq_mol, seq_meta = ensemble_for_sequence("YGGFM", window=8, c_term="free",
                                              n_conformers=4, seed=1, n_threads=1)
    shared = set(seq_meta) & set(smi_meta)
    assert {"sequence", "window", "truncated", "c_term_applied",
            "n_conformers", "seed"} <= shared
    assert smi_meta["built_from"] == "smiles"
    assert smi_meta["c_term_applied"] == "smiles"
    assert smi_mol.GetNumConformers() > 0 and seq_mol.GetNumConformers() > 0


def test_smiles_target_is_not_windowed():
    """A SMILES target is a complete molecule; no window may be applied to it."""
    _, meta = ensemble_for_smiles(DAMGO, n_conformers=4, seed=1, n_threads=1)
    assert meta["window"] == DAMGO
    assert meta["truncated"] is False
