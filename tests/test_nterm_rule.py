from rdkit import Chem

from mor_pharmacophore.features import peptide_mol
from mor_pharmacophore.nterm_rule import has_free_nterminus, supports_site


def test_free_peptide_has_free_nterminus():
    assert has_free_nterminus(peptide_mol("YGGFM"))


def test_blocked_nterminus_is_rejected():
    acetyl_ygg = Chem.MolFromSmiles("CC(=O)N[C@@H](Cc1ccc(O)cc1)C(=O)NCC(=O)O")
    assert not has_free_nterminus(acetyl_ygg)


def test_only_peptides_starting_at_the_tyr_support_the_site():
    assert supports_site(peptide_start=51, site_position=51)
    assert not supports_site(peptide_start=49, site_position=51)  # spans the Tyr
