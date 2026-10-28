#!/usr/bin/env python3
"""Is the screen deterministic? Separate the two places it could fail to be.

Written to settle a specific open question. When template A's dietary runs were
regenerated after the scratch allocation was reaped, 27 of 40 peptide-scores reproduced
to four decimal places but two did not: beta-lactorphin moved +0.0277 and casoxin A
+0.0213, roughly ten times the documented seed jitter, with conformer counts unchanged
(793->793, 800->800) so it is not an ensemble-size effect
(docs/dual_template_result.md, "A reproducibility note on A8").

Three candidate causes, and they need different fixes:

  C1  conformer generation is not deterministic given a fixed seed
  C2  align-it's alignment search is not deterministic given a fixed input
  C3  neither; the two runs used different library builds, because the environment was
      rebuilt from scratch between them

This script tests C1 and C2 directly, in one environment, by running each step twice on
identical input. C3 cannot be tested retrospectively -- the original run did not record
its library versions, which is itself the finding, and 04_screen.py now records them
(provenance block in `main`) so that this question is answerable next time.

A clean pass on C1 and C2 does not prove C3; it removes the alternatives, which is the
most this evidence can do. Say that in the write-up rather than claiming C3 confirmed.

Usage:
    python scripts/31_determinism_check.py --out docs/determinism_check.json
"""
import argparse
import csv
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from rdkit import Chem, RDLogger  # noqa: E402

from mor_pharmacophore.alignit import (  # noqa: E402
    babel_env, read_scores, score, version)
from mor_pharmacophore.conformers import ensemble_for_sequence  # noqa: E402

RDLogger.DisableLog("rdApp.*")

# The two peptides that failed to reproduce, plus one that did, as a negative control:
# a test that calls everything non-deterministic is not measuring anything.
#
# Keyed by `id` in the production targets table, with sequence AND c_term read from that
# table rather than restated here. The first version of this script hardcoded
# c_term="auto" and so built the wrong molecule for beta-lactorphin, which the table
# marks `amide`: 788 conformers against the 793 that both production runs produced. The
# determinism verdict was unaffected -- identical input twice is identical either way --
# but one of the two peptides under investigation was not the molecule the screen scored,
# which is the whole point of naming those two.
PEPTIDE_IDS = {
    "diet_beta_lactorphin": "moved +0.0277 on regeneration",
    "diet_casoxin_a": "moved +0.0213 on regeneration",
    "diet_met_enkephalin": "reproduced exactly; control",
}
DEFAULT_TARGETS = Path("data/dietary/dietary_targets.tsv")
SEED = 61453          # the seed whose run showed the discrepancy
WINDOW = 8
N_CONFORMERS = 800


def load_peptides(targets_path):
    """Sequence and c_term exactly as the production screen reads them."""
    rows = {}
    with Path(targets_path).open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["id"] in PEPTIDE_IDS:
                c_term = row.get("c_term") or "auto"
                if c_term not in ("auto", "free", "amide"):   # same guard as 04_screen
                    c_term = "auto"
                rows[row["id"]] = {"sequence": row["sequence"], "c_term": c_term,
                                   "why": PEPTIDE_IDS[row["id"]]}
    missing = set(PEPTIDE_IDS) - set(rows)
    if missing:
        raise SystemExit(f"targets table {targets_path} is missing {sorted(missing)}")
    return rows


def conformer_fingerprint(mol):
    """Hash every atom coordinate of every conformer, at full precision.

    Deliberately stricter than the 4-decimal comparison that surfaced the discrepancy:
    if embedding differs at all, this sees it.
    """
    h = hashlib.sha256()
    for conf in mol.GetConformers():
        pos = conf.GetPositions()
        h.update(pos.tobytes())
    return h.hexdigest(), mol.GetNumConformers()


def build_once(sequence, c_term):
    mol, info = ensemble_for_sequence(
        sequence, window=WINDOW, c_term=c_term,
        n_conformers=N_CONFORMERS, seed=SEED)
    digest, n = conformer_fingerprint(mol)
    return mol, info, digest, n


def write_sdf(mol, name, path):
    writer = Chem.SDWriter(str(path))
    mol.SetProp("_Name", name)
    for conf in mol.GetConformers():
        writer.write(mol, confId=conf.GetId())
    writer.close()
    return path


def protonate(in_sdf, out_sdf, ph=7.4):
    subprocess.run(["obabel", str(in_sdf), "-O", str(out_sdf), "-p", str(ph)],
                   env=babel_env(), capture_output=True, text=True, check=False)
    if not Path(out_sdf).exists() or Path(out_sdf).stat().st_size == 0:
        raise RuntimeError(f"obabel produced nothing for {in_sdf}")
    return out_sdf


def score_once(template, sdf, workdir, tag):
    out = Path(workdir) / f"scores_{tag}.tab"
    score(template, sdf, out, db_type="SDF", rank_by="TVERSKY_REF", with_exclusion=True)
    rows = read_scores(out)
    # Keep the full per-conformer vector: a mean can be stable while the alignment of
    # individual conformers is not, and that distinction is the point of C2.
    values = [r["tversky_ref"] for r in rows]
    return values


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--template",
                    default=Path("results/template_p8/mor_template_excl_a1.phar"),
                    type=Path)
    ap.add_argument("--targets", default=DEFAULT_TARGETS, type=Path,
                    help="production targets table; sequence and c_term are read "
                         "from it so the molecule matches what the screen scored")
    ap.add_argument("--out", default=Path("docs/determinism_check.json"), type=Path)
    args = ap.parse_args()

    if not args.template.exists():
        raise SystemExit(f"template not found: {args.template}")

    import rdkit
    env_record = {
        "rdkit": rdkit.__version__,
        "alignit": version(),
        "obabel": subprocess.run(["obabel", "-V"], env=babel_env(), capture_output=True,
                                 text=True).stdout.strip(),
        "python": sys.version.split()[0],
    }
    print("environment:", json.dumps(env_record, indent=2))

    peptides = load_peptides(args.targets)
    print("peptides (sequence and c_term read from the production targets table):",
          json.dumps(peptides, indent=2))

    results = {}
    workdir = Path(tempfile.mkdtemp(prefix="determinism_"))
    try:
        for name, spec in peptides.items():
            seq, c_term = spec["sequence"], spec["c_term"]
            print(f"\n=== {name} ({seq}, c_term={c_term}) -- {spec['why']} ===")

            # --- C1: conformer generation, twice, same seed --------------------------
            mol_a, info_a, dig_a, n_a = build_once(seq, c_term)
            mol_b, info_b, dig_b, n_b = build_once(seq, c_term)
            c1_identical = (dig_a == dig_b) and (n_a == n_b)
            print(f"  C1 conformers: n={n_a} vs {n_b}, "
                  f"{'IDENTICAL' if c1_identical else 'DIFFER'}")

            # --- C2: align-it, twice, on one fixed protonated SDF --------------------
            # One SDF scored twice, so any difference is the alignment search alone and
            # cannot be conformer variation leaking in.
            sdf = write_sdf(mol_a, name, workdir / f"{name}.sdf")
            prot = protonate(sdf, workdir / f"{name}_h.sdf")
            vals_a = score_once(args.template, prot, workdir, f"{name}_a")
            vals_b = score_once(args.template, prot, workdir, f"{name}_b")
            c2_identical = vals_a == vals_b
            max_delta = (max(abs(x - y) for x, y in zip(vals_a, vals_b))
                         if len(vals_a) == len(vals_b) and vals_a else None)
            print(f"  C2 align-it: {len(vals_a)} vs {len(vals_b)} conformers scored, "
                  f"{'IDENTICAL' if c2_identical else 'DIFFER'}"
                  + (f", max |delta| {max_delta:.6f}" if max_delta is not None else ""))

            mean_a = sum(vals_a) / len(vals_a) if vals_a else None
            mean_b = sum(vals_b) / len(vals_b) if vals_b else None
            results[name] = {
                "sequence": seq,
                "c_term": c_term,
                "why_included": spec["why"],
                "C1_conformers_identical": c1_identical,
                "C1_digest_run1": dig_a, "C1_digest_run2": dig_b,
                "C1_n_conformers": [n_a, n_b],
                "C2_alignit_identical": c2_identical,
                "C2_n_scored": [len(vals_a), len(vals_b)],
                "C2_max_abs_delta": max_delta,
                "tversky_ref_mean": [mean_a, mean_b],
            }
    finally:
        shutil.rmtree(workdir, ignore_errors=True)

    c1_all = all(r["C1_conformers_identical"] for r in results.values())
    c2_all = all(r["C2_alignit_identical"] for r in results.values())
    verdict = (
        "C1 and C2 both clean: within one environment the pipeline is bit-for-bit "
        "reproducible, so the regeneration discrepancy is not non-determinism. C3 (a "
        "different library build) is the remaining candidate and cannot be confirmed "
        "retrospectively, because the original run recorded no library versions."
        if c1_all and c2_all else
        f"NOT deterministic: C1 clean={c1_all}, C2 clean={c2_all}. The seed does not "
        f"pin the result, so every S1 seed-stability gate in this project is measuring "
        f"run-to-run noise as well as seed-to-seed noise."
    )
    print("\n" + verdict)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({
        "generated": str(date.today()),
        "question": "is the screen deterministic given a fixed seed and fixed input?",
        "seed": SEED, "window": WINDOW, "n_conformers": N_CONFORMERS,
        "template": str(args.template),
        "environment": env_record,
        "peptides": results,
        "C1_conformer_generation_deterministic": c1_all,
        "C2_alignit_deterministic": c2_all,
        "verdict": verdict,
    }, indent=2) + "\n")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
