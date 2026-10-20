#!/usr/bin/env python3
"""Assemble the single list of everything the screen scores.

The library and the controls are scored in one pass, under identical settings. Scoring
controls in a separate run would let a difference in protonation, conformer count or
Align-it invocation masquerade as a difference in binding.

Each row carries a `set` label so the validation stage can partition the scores without
re-deriving membership:

    library             a digestion peptide
    positive_endogenous an IUPHAR/BPS endogenous ligand of OPRM1
    positive_food       a documented food-derived MOR ligand
    template_ligand     met-enkephalin, the ligand the template was read from

Set membership is not exclusive in principle - beta-casomorphin-7 is both a library
member and a documented active - so the label records why a peptide is being scored, and
the library row is kept as the canonical one.

Usage:
    python scripts/03b_targets.py --out data/targets.tsv
"""
import argparse
import csv
import json
from pathlib import Path


def read_tsv(path):
    with Path(path).open() as fh:
        rows = [l for l in fh if not l.startswith("#") and l.strip()]
    return list(csv.DictReader(rows, delimiter="\t"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--library", default=Path("data/library/peptides.tsv"), type=Path)
    ap.add_argument("--endogenous",
                    default=Path("data/reference/oprm1_endogenous_ligands.tsv"), type=Path)
    ap.add_argument("--reference",
                    default=Path("data/reference/reference_peptides.tsv"), type=Path)
    ap.add_argument("--out", default=Path("data/targets.tsv"), type=Path)
    args = ap.parse_args()

    targets, seen = [], set()

    for row in read_tsv(args.library):
        targets.append({"id": row["peptide_id"], "sequence": row["sequence"],
                        "set": "library", "c_term": "auto"})
        seen.add(row["peptide_id"])

    for row in read_tsv(args.endogenous):
        label = "template_ligand" if row["use"] == "template" else "positive_endogenous"
        tid = f"ref__{row['name']}"
        if tid in seen:
            continue
        targets.append({"id": tid, "sequence": row["sequence"], "set": label,
                        "c_term": row["c_term"]})
        seen.add(tid)

    skipped = []
    for row in read_tsv(args.reference):
        if row["sequence"] == "NON_CANONICAL":
            skipped.append(row["name"])       # synthetic, no plain sequence to model
            continue
        tid = f"ref__{row['name']}"
        if tid in seen:
            continue
        c_term = row["c_term"] if row["c_term"] in ("free", "amide") else "auto"
        targets.append({"id": tid, "sequence": row["sequence"],
                        "set": "positive_food", "c_term": c_term})
        seen.add(tid)

    columns = ["id", "sequence", "set", "c_term"]
    with args.out.open("w") as fh:
        fh.write("\t".join(columns) + "\n")
        for t in targets:
            fh.write("\t".join(t[c] for c in columns) + "\n")

    counts = {}
    for t in targets:
        counts[t["set"]] = counts.get(t["set"], 0) + 1
    summary = {"total": len(targets), "by_set": counts, "skipped_non_canonical": skipped}
    args.out.with_suffix(".summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"wrote {args.out}")
    for k, v in counts.items():
        print(f"  {k:24s} {v}")
    if skipped:
        print(f"  skipped (non-canonical): {skipped}")


if __name__ == "__main__":
    main()
