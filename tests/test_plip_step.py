from pathlib import Path

import pytest

from mor_pharmacophore.plip_step import peptide_interactions

COMPLEX = Path(__file__).parent.parent / "data" / "8f7q_complex.pdb"


@pytest.fixture(scope="module")
def ix():
    return peptide_interactions(COMPLEX, "P")


def test_peptide_recovered_as_ligand(ix):
    # chain P of 8F7Q residues 1-5 is the Tyr-Gly-Gly-Phe-Met (met-enkephalin) segment
    assert ix.ligand == "TYR-GLY-GLY-PHE-MET"
    assert ix.chain == "P" and ix.position == 1


def test_interactions_come_from_receptor_only(ix):
    residues = (ix.hbond_residues + ix.hydrophobic_residues
                + ix.saltbridge_residues + ix.pistack_residues)
    assert residues and all(r.endswith(":R") for r in residues)


def test_hydrogen_bonds_and_hydrophobic_pocket(ix):
    # regression snapshot of the PLIP result on the experimental pose
    assert len(ix.hbond_residues) >= 2
    assert len(ix.hydrophobic_residues) >= 4
    assert {"TYR150:R", "TRP320:R"} <= set(ix.hydrophobic_residues)
