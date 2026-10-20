#!/usr/bin/env python3
"""Stage 3b: freeze the control sets, before any screening run.

Four sets, all derived here so that none can be chosen after seeing a score
(PLAN.md 5):

  positive_endogenous  IUPHAR/BPS-curated endogenous ligands of OPRM1.
  positive_food        documented food-derived MOR ligands (curated reference set).
  positional_control   library peptides containing the YPFP message, split by whether
                       the peptide STARTS with it (free N-terminal Tyr) or merely SPANS
                       it (Tyr amine locked in a peptide bond). Same motif, same parent
                       protein, different position - the set that isolates the
                       positional rule from motif presence.
  baseline_matched     Tyr-starting library peptides that are not known actives. The
                       sequence baseline "starts with Tyr" is constant across this set,
                       so any separation within it is attributable to 3D matching.

Usage:
    python scripts/02_controls.py --library data/library/peptides.tsv --out data/controls
"""
import argparse
import csv
import json
from datetime import date
from pathlib import Path

# The message motif of the beta-casomorphins. Chosen because the library contains both
# peptides that start with it and many more that span it, which is what makes the
# positional control possible at all.
POSITIONAL_MOTIF = "YPFP"

REFERENCE = Path("data/reference/reference_peptides.tsv")
ENDOGENOUS = Path("data/reference/oprm1_endogenous_ligands.tsv")


def read_tsv(path):
    with path.open() as fh:
        rows = [l for l in fh if not l.startswith("#") and l.strip()]
    return list(csv.DictReader(rows, delimiter="\t"))


def read_library(path):
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--library", default=Path("data/library/peptides.tsv"), type=Path)
    ap.add_argument("--out", default=Path("data/controls"), type=Path)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    library = read_library(args.library)
    by_sequence = {r["sequence"]: r for r in library}

    endogenous = [r for r in read_tsv(ENDOGENOUS) if r["use"] == "positive"]
    template_ligand = [r for r in read_tsv(ENDOGENOUS) if r["use"] == "template"]

    # Food-derived actives. CTOP is a non-canonical synthetic peptide with no plain
    # sequence, so it cannot enter a sequence-driven screen and is dropped by name.
    food = [r for r in read_tsv(REFERENCE) if r["sequence"] != "NON_CANONICAL"]

    sets = {}

    sets["positive_endogenous"] = [
        {"name": r["name"], "sequence": r["sequence"], "c_term": r["c_term"],
         "role": r["role"], "pKi": r["pKi"], "pmid": r["pmid"],
         "in_library": r["sequence"] in by_sequence}
        for r in endogenous
    ]

    sets["positive_food"] = [
        {"name": r["name"], "sequence": r["sequence"], "c_term": r["c_term"],
         "role": r["role"], "starts_with_Y": r["sequence"].startswith("Y"),
         "in_library": r["sequence"] in by_sequence}
        for r in food
    ]

    starts, spans = [], []
    for r in library:
        seq = r["sequence"]
        if POSITIONAL_MOTIF not in seq:
            continue
        entry = {"peptide_id": r["peptide_id"], "sequence": seq,
                 "length": int(r["length"]), "phases": r["phases"]}
        (starts if seq.startswith(POSITIONAL_MOTIF) else spans).append(entry)
    sets["positional_control"] = {
        "motif": POSITIONAL_MOTIF,
        "starts_with_motif": starts,
        "spans_motif_blocked_nterm": spans,
    }

    active_sequences = {r["sequence"] for r in endogenous} | {r["sequence"] for r in food}
    active_sequences |= {r["sequence"] for r in template_ligand}
    sets["baseline_matched"] = [
        {"peptide_id": r["peptide_id"], "sequence": r["sequence"],
         "length": int(r["length"]), "phases": r["phases"]}
        for r in library
        if r["starts_with_Y"] == "yes" and r["sequence"] not in active_sequences
    ]

    sets["template_ligand"] = template_ligand

    (args.out / "control_sets.json").write_text(json.dumps(sets, indent=2) + "\n")

    summary = {
        "generated": date.today().isoformat(),
        "library": str(args.library),
        "positive_endogenous": len(sets["positive_endogenous"]),
        "positive_endogenous_in_library": sum(
            1 for r in sets["positive_endogenous"] if r["in_library"]),
        "positive_food": len(sets["positive_food"]),
        "positive_food_in_library": sum(
            1 for r in sets["positive_food"] if r["in_library"]),
        "positive_food_not_Y_start": [
            r["name"] for r in sets["positive_food"] if not r["starts_with_Y"]],
        "positional_starts": len(starts),
        "positional_spans": len(spans),
        "baseline_matched": len(sets["baseline_matched"]),
    }
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    for k, v in summary.items():
        print(f"  {k:34s} {v}")


if __name__ == "__main__":
    main()
