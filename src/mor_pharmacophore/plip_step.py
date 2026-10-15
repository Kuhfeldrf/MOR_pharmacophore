"""Step 1: read the interactions a receptor makes with its own peptide ligand (PLIP)."""
from collections import defaultdict
from dataclasses import dataclass, field

from plip.basic import config
from plip.structure.preparation import PDBComplex


@dataclass
class PeptideInteractions:
    ligand: str
    chain: str
    position: int
    hbond_residues: list = field(default_factory=list)
    hydrophobic_residues: list = field(default_factory=list)
    saltbridge_residues: list = field(default_factory=list)
    pistack_residues: list = field(default_factory=list)

    @property
    def n_interactions(self):
        return (len(self.hbond_residues) + len(self.hydrophobic_residues)
                + len(self.saltbridge_residues) + len(self.pistack_residues))


def _res(i):
    return f"{i.restype}{i.resnr}:{i.reschain}"


def peptide_interactions(pdb_path, peptide_chain):
    """Run PLIP with `peptide_chain` treated as the ligand; return its interactions."""
    config.PEPTIDES = [peptide_chain]  # PLIP ignores peptide ligands unless told which chain
    mol = PDBComplex()
    mol.load_pdb(str(pdb_path))
    mol.analyze()
    found = [lig for lig in mol.ligands if lig.chain == peptide_chain]
    if not found:
        raise ValueError(f"PLIP found no peptide ligand on chain {peptide_chain}")
    lig = found[0]
    s = mol.interaction_sets[":".join([lig.hetid, lig.chain, str(lig.position)])]
    return PeptideInteractions(
        ligand=lig.longname, chain=lig.chain, position=lig.position,
        hbond_residues=sorted({_res(i) for i in list(s.hbonds_pdon) + list(s.hbonds_ldon)}),
        hydrophobic_residues=sorted({_res(i) for i in s.hydrophobic_contacts}),
        saltbridge_residues=sorted({_res(i) for i in list(s.saltbridge_lneg) + list(s.saltbridge_pneg)}),
        pistack_residues=sorted({_res(i) for i in s.pistacking}),
    )
