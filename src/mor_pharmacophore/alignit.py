"""Step 3: Align-it (external binary, GPL/LGPL-3.0, https://github.com/OliverBScott/align-it).

Align-it both generates a pharmacophore from a molecule and scores candidates against
one. It is not on conda-forge and must be built from source; see
`env/build_alignit.sbatch`. The build needs the conda toolchain rather than a system
g++, because conda-forge OpenBabel requires a newer CXXABI than RHEL's libstdc++ carries.

Feature codes used by Align-it: AROM, HDON, HACC, LIPO, CHARGE.
"""
import shutil
import subprocess
from pathlib import Path

# Upstream installs the binary as `align_it`; some packagings name it `align-it`.
_CANDIDATE_NAMES = ("align_it", "align-it")

FEATURE_CODES = ("AROM", "HDON", "HACC", "LIPO", "CHARGE")


def alignit_path():
    """Absolute path to the Align-it binary, or None if it is not on PATH."""
    for name in _CANDIDATE_NAMES:
        found = shutil.which(name)
        if found:
            return found
    return None


def _require():
    exe = alignit_path()
    if exe is None:
        raise FileNotFoundError(
            "align-it is not on PATH (looked for: " + ", ".join(_CANDIDATE_NAMES)
            + "). Build it with env/build_alignit.sbatch.")
    return exe


def _run(args):
    result = subprocess.run([_require()] + args, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"align-it failed ({result.returncode}): {result.stderr.strip()[:500]}")
    return result


def alignit_help():
    return subprocess.run([_require(), "--help"], capture_output=True, text=True)


def version():
    """Version banner. Align-it prints its banner on stderr, not stdout."""
    result = alignit_help()
    for line in (result.stderr + result.stdout).splitlines():
        if "Align-it v" in line:
            return line.split("|")[0].replace("Align-it", "align-it").strip()
    return "unknown"


def generate_pharmacophore(molecule_path, out_phar, mol_type=None):
    """Generate the pharmacophore of `molecule_path` and write it to `out_phar`.

    This is how the template is built from the experimental ligand pose: Align-it reads
    the 3D structure and emits typed feature points at their coordinates.
    """
    args = ["-d", str(molecule_path), "-p", str(out_phar), "--noHybrid"]
    if mol_type:
        args += ["--dbType", mol_type]
    _run(args)
    return Path(out_phar)


def score(reference_phar, dbase_path, scores_out, aligned_out=None,
          rank_by="TANIMOTO", cutoff=None, best=None, db_type=None):
    """Score every molecule in `dbase_path` against the `reference_phar` template.

    Returns the path to the score table. `rank_by` TVERSKY_REF weights the reference's
    own features, which is the asymmetric option to use when candidates are larger than
    the template.
    """
    args = ["-r", str(reference_phar), "--refType", "PHAR",
            "-d", str(dbase_path), "-s", str(scores_out), "--rankBy", rank_by]
    if aligned_out:
        args += ["-o", str(aligned_out)]
    if db_type:
        args += ["--dbType", db_type]
    if cutoff is not None:
        args += ["--cutOff", str(cutoff)]
    if best is not None:
        args += ["--best", str(best)]
    _run(args)
    return Path(scores_out)


# Column layout of the Align-it score table, verified empirically against v1.0.4 output
# (headerless, tab-separated, one row per database conformer):
#   1 reference name      2 reference volume   3 candidate name   4 candidate volume
#   5 overlap             6 exclusion-volume overlap              7 corrected overlap
#   8 matched features    9 TANIMOTO          10 TVERSKY_REF     11 TVERSKY_DB
_SCORE_COLUMNS = 11


def read_scores(scores_path):
    """Parse an Align-it score table into dicts, one per database conformer."""
    rows = []
    for line in Path(scores_path).read_text().splitlines():
        if not line.strip():
            continue
        parts = line.split("\t") if "\t" in line else line.split()
        if len(parts) < _SCORE_COLUMNS:
            continue
        rows.append({
            "reference": parts[0],
            "reference_volume": float(parts[1]),
            "candidate": parts[2],
            "candidate_volume": float(parts[3]),
            "overlap": float(parts[4]),
            "corrected_overlap": float(parts[6]),
            "matched_features": int(float(parts[7])),
            "tanimoto": float(parts[8]),
            "tversky_ref": float(parts[9]),
            "tversky_db": float(parts[10]),
        })
    return rows


def best_per_candidate(rows, metric="tversky_ref"):
    """Collapse a conformer-level score table to one row per candidate.

    The per-peptide score is the best over its ensemble; the spread across the ensemble
    is retained, because cross-conformer agreement is the convergence readout (C5) and
    in this project's earlier work it was the only signal that separated a good pose
    from a bad one.
    """
    grouped = {}
    for row in rows:
        grouped.setdefault(row["candidate"], []).append(row)
    out = {}
    for name, group in grouped.items():
        values = sorted(r[metric] for r in group)
        best = max(group, key=lambda r: r[metric])
        n = len(values)
        mean = sum(values) / n
        out[name] = {
            "candidate": name,
            "n_conformers": n,
            metric: best[metric],
            "tanimoto": best["tanimoto"],
            "tversky_ref": best["tversky_ref"],
            "tversky_db": best["tversky_db"],
            "matched_features": best["matched_features"],
            f"{metric}_mean": mean,
            f"{metric}_min": values[0],
            f"{metric}_spread": values[-1] - values[0],
        }
    return out
