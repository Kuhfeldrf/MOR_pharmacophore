"""Stage 2: build 3D conformer ensembles for candidate peptides.

Two decisions govern this module, both recorded in PLAN.md and both load-bearing.

**D3 - the N-terminal window.** Only the first `WINDOW` residues are modelled. The
library's median length is 15 and its maximum 46; the template is 5 residues. Extra
residues contribute pharmacophore points with nothing to match, which dilutes the score
- the same failure already measured in this project as ipTM dilution, where 281 residues
against 5 was called "close to the worst case". The window is also what the binding
model describes: the opioid message is an N-terminal motif anchored by the free
alpha-amine.

**D5 - the truncation cap.** When a peptide is longer than the window, the window's
C-terminus is an artefact of cutting, not a real chain end. Leaving it as a free
carboxylate would add a negative charge and acceptor features that the intact peptide
does not have at that position. Truncated windows are therefore C-terminally amidated.
A peptide that genuinely ends within the window keeps its real C-terminus (free acid, or
amide where the reference data records one).
"""
from rdkit import Chem
from rdkit.Chem import AllChem

WINDOW = 5
N_CONFORMERS = 50
PRUNE_RMS = 0.5
MMFF_MAX_ITERS = 500


def _amidate_c_terminus(mol):
    """Replace the C-terminal carboxylate/carboxylic acid with a primary amide."""
    patt = Chem.MolFromSmarts("[CX3](=O)[OX2H1,OX1-]")
    matches = mol.GetSubstructMatches(patt)
    if not matches:
        return mol
    # The C-terminus is the last such group in the molecule's atom ordering.
    _, _, o_idx = matches[-1]
    edit = Chem.RWMol(mol)
    atom = edit.GetAtomWithIdx(o_idx)
    atom.SetAtomicNum(7)
    atom.SetFormalCharge(0)
    atom.SetNoImplicit(False)
    atom.SetNumExplicitHs(2)
    out = edit.GetMol()
    Chem.SanitizeMol(out)
    return out


def build_window(sequence, window=WINDOW, c_term="auto"):
    """Return (mol, was_truncated) for the N-terminal window of `sequence`.

    c_term: "auto"   - amidate if truncated, otherwise free acid (D5)
            "amide"  - amidate regardless (peptide is genuinely amidated)
            "free"   - leave the free acid
    """
    if not sequence:
        raise ValueError("empty sequence")
    truncated = len(sequence) > window
    fragment = sequence[:window]

    mol = Chem.MolFromSequence(fragment)
    if mol is None:
        raise ValueError(f"RDKit could not build a peptide from {fragment!r}")

    if c_term == "amide" or (c_term == "auto" and truncated):
        mol = _amidate_c_terminus(mol)

    return mol, truncated


def embed_conformers(mol, n_conformers=N_CONFORMERS, seed=0xF00D,
                     prune_rms=PRUNE_RMS, minimise=True):
    """Embed and MMFF-minimise a conformer ensemble. Returns the mol with conformers."""
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.pruneRmsThresh = prune_rms
    params.useSmallRingTorsions = True
    ids = AllChem.EmbedMultipleConfs(mol, numConfs=n_conformers, params=params)
    if len(ids) == 0:
        raise ValueError("no conformers could be embedded")
    if minimise:
        # Failures are tolerated per-conformer: MMFF has no parameters for some states,
        # and dropping the ensemble over one conformer would lose the peptide.
        try:
            AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=MMFF_MAX_ITERS)
        except Exception:
            pass
    return mol


def ensemble_for_sequence(sequence, window=WINDOW, c_term="auto",
                          n_conformers=N_CONFORMERS, seed=0xF00D):
    """Convenience: window + embed in one call. Returns (mol, metadata)."""
    mol, truncated = build_window(sequence, window=window, c_term=c_term)
    mol = embed_conformers(mol, n_conformers=n_conformers, seed=seed)
    meta = {
        "sequence": sequence,
        "window": sequence[:window],
        "truncated": truncated,
        "c_term_applied": "amide" if (c_term == "amide" or (c_term == "auto" and truncated)) else "free",
        "n_conformers": mol.GetNumConformers(),
        "seed": seed,
    }
    return mol, meta


def write_sdf(mol, path, name=None, properties=None):
    """Write every conformer of `mol` as a record in one SDF."""
    if name:
        mol.SetProp("_Name", name)
    if properties:
        for k, v in properties.items():
            mol.SetProp(str(k), str(v))
    writer = Chem.SDWriter(str(path))
    for conf in mol.GetConformers():
        writer.write(mol, confId=conf.GetId())
    writer.close()
    return path
