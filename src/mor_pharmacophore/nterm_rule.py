"""The positional requirement: opioid binding needs a FREE N-terminus at the anchoring Tyr."""
from rdkit import Chem

_FREE_NTERM = Chem.MolFromSmarts("[NX3;H2,H3+][CX4][CX3](=O)")  # primary amine on alpha carbon


def has_free_nterminus(mol):
    """True if the first residue's alpha-amine is unblocked (not acetylated, not an amide)."""
    return mol.HasSubstructMatch(_FREE_NTERM)


def supports_site(peptide_start, site_position):
    """A digestion peptide supports a Tyr site only if its N-terminus sits exactly there.

    Peptides that merely span the Tyr do not count: their N-terminus is further upstream
    and the Tyr amine is blocked in the peptide bond.
    """
    return peptide_start == site_position
