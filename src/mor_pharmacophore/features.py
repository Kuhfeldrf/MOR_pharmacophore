"""Step 2: encode a peptide as RDKit pharmacophore features (BaseFeatures.fdef)."""
import os
from collections import Counter

from rdkit import Chem, RDConfig
from rdkit.Chem import AllChem, ChemicalFeatures

_FACTORY = ChemicalFeatures.BuildFeatureFactory(
    os.path.join(RDConfig.RDDataDir, "BaseFeatures.fdef"))


def peptide_mol(sequence, seed=7):
    """3D model of a free (unblocked) peptide from one-letter sequence."""
    mol = Chem.AddHs(Chem.MolFromSequence(sequence))
    if AllChem.EmbedMolecule(mol, randomSeed=seed) != 0:
        raise ValueError(f"could not embed {sequence}")
    return mol


def feature_counts(sequence):
    """Count pharmacophore features by family, e.g. {'Aromatic': 2, 'Donor': 7}."""
    mol = peptide_mol(sequence)
    return Counter(f.GetFamily() for f in _FACTORY.GetFeaturesForMol(mol))


def feature_positions(sequence, family):
    """3D positions of one feature family (the template Align-it superimposes onto)."""
    mol = peptide_mol(sequence)
    out = []
    for f in _FACTORY.GetFeaturesForMol(mol):
        if f.GetFamily() == family:
            p = f.GetPos()
            out.append((p.x, p.y, p.z))
    return out
