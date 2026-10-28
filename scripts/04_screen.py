#!/usr/bin/env python3
"""Stages 2+4: conformers and Align-it scoring for one block of peptides.

One block = one array task. Blocks are independent, so the array is embarrassingly
parallel and the number running at once is set by the scheduler, not by this script.

Within a block every peptide's conformers go into a single SDF, protonated in one
OpenBabel call and scored in one Align-it call. Protonation is applied the same way as
for the template (pH 7.4), so template and candidates are on the same footing - without
that the N-terminal ammonium anchor would be present on one side only.

Usage:
    python scripts/04_screen.py --targets data/targets.tsv --template <phar> \\
        --block 0 --nblocks 146 --out results/04_screen
"""
import argparse
import csv
import json
import sys
import subprocess
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from rdkit import Chem, RDLogger  # noqa: E402

from mor_pharmacophore.alignit import (  # noqa: E402
    babel_env, best_per_candidate, read_scores, score, version)
from mor_pharmacophore.conformers import (  # noqa: E402
    N_CONFORMERS, WINDOW, ensemble_for_sequence, ensemble_for_smiles)

RDLogger.DisableLog("rdApp.*")


def environment_record():
    """The library stack that produced these scores.

    Deliberately tolerant: a missing version must degrade this record, never fail a
    screen that has already spent its CPU. Unknown entries are recorded as null so the
    gap is visible rather than absent.
    """
    import platform

    import rdkit

    def safe(fn):
        """Every lookup is independently fallible; one missing tool loses one field."""
        try:
            return fn() or None
        except Exception:
            return None

    return {
        "rdkit": getattr(rdkit, "__version__", None),
        "alignit": safe(version),
        "obabel": safe(lambda: subprocess.run(
            ["obabel", "-V"], env=babel_env(),
            capture_output=True, text=True).stdout.strip()),
        "python": platform.python_version(),
        "platform": platform.platform(),
    }


def read_targets(path):
    with Path(path).open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def block_slice(rows, block, nblocks):
    """Deterministic contiguous slice, so a rerun of block k covers the same peptides."""
    if nblocks <= 1:
        return rows
    size = (len(rows) + nblocks - 1) // nblocks
    return rows[block * size:(block + 1) * size]


def build_block_sdf(rows, sdf_path, n_conformers, window, seed):
    """Embed every peptide in the block into one SDF. Returns (metadata, failures)."""
    writer = Chem.SDWriter(str(sdf_path))
    meta, failures = {}, []
    for row in rows:
        name, sequence = row["id"], row["sequence"]
        c_term = row.get("c_term") or "auto"
        if c_term not in ("auto", "free", "amide"):
            c_term = "auto"
        # An explicit `smiles` overrides `sequence`, for ligands whose backbone chemistry
        # has no one-letter code (DAMGO's N-MePhe and glycinol). Same embedding protocol
        # either way, so the scores stay comparable; see conformers.build_from_smiles.
        smiles = (row.get("smiles") or "").strip()
        try:
            if smiles:
                mol, info = ensemble_for_smiles(
                    smiles, n_conformers=n_conformers, seed=seed)
            else:
                mol, info = ensemble_for_sequence(
                    sequence, window=window, c_term=c_term,
                    n_conformers=n_conformers, seed=seed)
        except Exception as exc:                      # recorded, never silently dropped
            failures.append({"id": name, "sequence": sequence, "error": str(exc)[:200]})
            continue
        mol.SetProp("_Name", name)
        for conf in mol.GetConformers():
            writer.write(mol, confId=conf.GetId())
        meta[name] = info
    writer.close()
    return meta, failures


def protonate_sdf(in_sdf, out_sdf, ph=7.4):
    result = subprocess.run(
        ["obabel", str(in_sdf), "-O", str(out_sdf), "-p", str(ph)],
        env=babel_env(),   # pin OpenBabel's plugin dir; see alignit.babel_env
        capture_output=True, text=True)
    if not Path(out_sdf).exists() or Path(out_sdf).stat().st_size == 0:
        raise RuntimeError(f"obabel protonation produced nothing: {result.stderr[:300]}")
    return out_sdf


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", default=Path("data/targets.tsv"), type=Path)
    ap.add_argument("--template", default=Path("results/01_template/mor_template.phar"),
                    type=Path)
    ap.add_argument("--out", default=Path("results/04_screen"), type=Path)
    ap.add_argument("--block", default=0, type=int)
    ap.add_argument("--nblocks", default=1, type=int)
    ap.add_argument("--n-conformers", default=N_CONFORMERS, type=int)
    ap.add_argument("--window", default=WINDOW, type=int)
    ap.add_argument("--seed", default=0xF00D, type=int)
    ap.add_argument("--with-exclusion", action="store_true",
                    help="honour EXCL spheres in the template (align-it ignores them "
                         "otherwise, so a template with exclusion volumes scored "
                         "without this gives identical numbers to one without)")
    ap.add_argument("--keep-sdf", action="store_true",
                    help="retain conformer SDFs (they are transient by default; "
                         "scratch is tight)")
    args = ap.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    tag = f"block{args.block:04d}"
    work = args.out / tag
    work.mkdir(exist_ok=True)

    rows = block_slice(read_targets(args.targets), args.block, args.nblocks)
    print(f"{tag}: {len(rows)} peptides, align-it {version()}, "
          f"window {args.window}, {args.n_conformers} conformers, seed {args.seed}")

    t0 = time.time()
    raw_sdf = work / "conformers.sdf"
    meta, failures = build_block_sdf(
        rows, raw_sdf, args.n_conformers, args.window, args.seed)
    t_conf = time.time() - t0
    print(f"  conformers: {len(meta)} ok, {len(failures)} failed, {t_conf:.1f}s")

    if not meta:
        raise SystemExit(f"{tag}: no conformers built; refusing to write an empty score "
                         f"table ({len(failures)} failures)")

    t1 = time.time()
    prot_sdf = protonate_sdf(raw_sdf, work / "conformers_h.sdf")
    scores_tab = work / "scores.tab"
    score(args.template, prot_sdf, scores_tab, db_type="SDF", rank_by="TVERSKY_REF",
          with_exclusion=args.with_exclusion)
    t_score = time.time() - t1
    print(f"  scoring: {t_score:.1f}s")

    best = best_per_candidate(read_scores(scores_tab), metric="tversky_ref")

    by_id = {r["id"]: r for r in rows}
    out_tsv = args.out / f"{tag}.tsv"
    columns = ["id", "sequence", "set", "window", "truncated", "c_term_applied",
               "n_conformers", "tversky_ref", "tversky_ref_mean", "tversky_ref_spread",
               "tanimoto", "tversky_db", "matched_features"]
    with out_tsv.open("w") as fh:
        fh.write("\t".join(columns) + "\n")
        for name, info in sorted(best.items(), key=lambda kv: -kv[1]["tversky_ref"]):
            m = meta.get(name, {})
            src = by_id.get(name, {})
            fh.write("\t".join(str(v) for v in [
                name, src.get("sequence", ""), src.get("set", ""),
                m.get("window", ""), m.get("truncated", ""), m.get("c_term_applied", ""),
                info["n_conformers"], f"{info['tversky_ref']:.4f}",
                f"{info['tversky_ref_mean']:.4f}", f"{info['tversky_ref_spread']:.4f}",
                f"{info['tanimoto']:.4f}", f"{info['tversky_db']:.4f}",
                info["matched_features"],
            ]) + "\n")

    (args.out / f"{tag}.json").write_text(json.dumps({
        "block": args.block, "nblocks": args.nblocks,
        "n_targets": len(rows), "n_scored": len(best),
        "failures": failures,
        "alignit_version": version(),
        # Library versions, not just align-it's. Template A's dietary runs were
        # regenerated after the scratch allocation was reaped, two of twenty peptides
        # moved ~10x seed jitter, and the obvious explanation - the environment was
        # rebuilt from scratch in between - could not be checked, because the original
        # run recorded only align-it's version. A seed is not provenance if the stack
        # that consumed it is unrecorded. See docs/determinism_check.md.
        "environment": environment_record(),
        "window": args.window, "n_conformers": args.n_conformers, "seed": args.seed,
        "seconds_conformers": round(t_conf, 1), "seconds_scoring": round(t_score, 1),
    }, indent=2) + "\n")

    # The raw per-conformer table is retained: ensemble statistics (max / mean /
    # top-k) are read from it post hoc, so changing the statistic costs no recompute.
    if not args.keep_sdf:
        for f in (raw_sdf, prot_sdf):
            Path(f).unlink(missing_ok=True)

    print(f"  wrote {out_tsv} ({len(best)} scored)")


if __name__ == "__main__":
    main()
