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
import os

from rdkit import Chem
from rdkit.Chem import AllChem

def allocated_threads():
    """Threads this process may actually use.

    RDKit's `numThreads=0` means "every core the machine has", NOT "every core this job
    was allocated". Under a scheduler that packs several tasks onto a node, 0 therefore
    oversubscribes badly: eight 8-core tasks on a 64-core node each spawned 64 threads
    and drove the node to CPULoad 319 (measured 2026-10-20, sweep 237134). Default to
    the allocation instead, and let an explicit argument override it.
    """
    for var in ("SLURM_CPUS_PER_TASK", "OMP_NUM_THREADS"):
        value = os.environ.get(var)
        if value and value.isdigit() and int(value) > 0:
            return int(value)
    return len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else 1


WINDOW = 5
N_CONFORMERS = 50
PRUNE_RMS = 0.5
# Fraction of the requested conformers that must actually embed before the ensemble is
# accepted. Guards against the silent under-delivery documented in embed_conformers.
MIN_CONFORMER_FRACTION = 0.5
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
    # window=0 (or >= the length) means the whole peptide: no truncation, no cap.
    if window and window > 0:
        truncated = len(sequence) > window
        fragment = sequence[:window]
    else:
        truncated = False
        fragment = sequence

    mol = Chem.MolFromSequence(fragment)
    if mol is None:
        raise ValueError(f"RDKit could not build a peptide from {fragment!r}")

    if c_term == "amide" or (c_term == "auto" and truncated):
        mol = _amidate_c_terminus(mol)

    return mol, truncated


def embed_conformers(mol, n_conformers=N_CONFORMERS, seed=0xF00D,
                     prune_rms=PRUNE_RMS, minimise=True, n_threads=None,
                     use_random_coords=True, min_fraction=MIN_CONFORMER_FRACTION):
    """Embed and MMFF-minimise a conformer ensemble. Returns the mol with conformers.

    `n_threads=None` uses the job's allocation (see `allocated_threads`). Both the
    embedding and the minimisation are threaded; leaving them single-threaded is the
    dominant cost of the whole screen, and this project has already lost a campaign to
    an unthreaded run (OMP_NUM_THREADS=1, a measured 4.3x). Passing 0 here would mean
    "all hardware cores" and oversubscribes a shared node, so it is not the default.
    """
    if n_threads is None:
        n_threads = allocated_threads()
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.pruneRmsThresh = prune_rms
    params.useSmallRingTorsions = True
    params.numThreads = n_threads
    # useRandomCoords is REQUIRED, not an optimisation. Distance-geometry embedding
    # silently under-delivers on flexible peptides: measured 2026-10-20 on the
    # N-terminal window of YGGFLRRIRPKLKWDNQ, asking for 50 conformers returned
    #   window  5 (19 rotatable bonds) -> 48
    #   window  8 (37)                 -> 21
    #   window 10 (44)                 ->  2
    #   window 12 (57)                 ->  1
    # with no error raised. With random starting coordinates all of these return 50.
    # The first screen ran at window 5 and so never saw this; any longer window would
    # have scored tiny ensembles as though they were full ones.
    params.useRandomCoords = use_random_coords
    ids = AllChem.EmbedMultipleConfs(mol, numConfs=n_conformers, params=params)
    if len(ids) == 0:
        raise ValueError("no conformers could be embedded")
    # A short ensemble is a quiet way to get a meaningless best-of-ensemble score, so
    # it is surfaced rather than tolerated.
    if len(ids) < min_fraction * n_conformers:
        raise ValueError(
            f"only {len(ids)} of {n_conformers} conformers embedded "
            f"({len(ids) / n_conformers:.0%} < {min_fraction:.0%} required)")
    if minimise:
        # Failures are tolerated per-conformer: MMFF has no parameters for some states,
        # and dropping the ensemble over one conformer would lose the peptide.
        try:
            AllChem.MMFFOptimizeMoleculeConfs(
                mol, maxIters=MMFF_MAX_ITERS, numThreads=n_threads)
        except Exception:
            pass
    return mol


def ensemble_for_sequence(sequence, window=WINDOW, c_term="auto",
                          n_conformers=N_CONFORMERS, seed=0xF00D, n_threads=None,
                          use_random_coords=True):
    """Convenience: window + embed in one call. Returns (mol, metadata)."""
    mol, truncated = build_window(sequence, window=window, c_term=c_term)
    mol = embed_conformers(mol, n_conformers=n_conformers, seed=seed,
                           n_threads=n_threads,
                           use_random_coords=use_random_coords)
    meta = {
        "sequence": sequence,
        "window": sequence[:window] if (window and window > 0) else sequence,
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
