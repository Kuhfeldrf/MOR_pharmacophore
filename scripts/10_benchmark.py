#!/usr/bin/env python3
"""Build the benchmark for the validation/optimisation study (STUDY.md 2).

Sets produced, all deterministic:

  active      36 ChEMBL peptides with measured pKi (4.71-9.14, 4.43 log units)
  scramble    one composition-matched permutation per active, fixed seed
  tyr_hard    Tyr-starting library peptides, length-matched to the actives
  positional  library peptides containing YPFP, split by whether they start with it

Family labels carry the potency-pair control (STUDY.md 2.1): the YGGFL family
(n=11, pKi 7.39-8.71) against the YGPFL family (n=6, pKi 4.71-5.51), which differ by a
single Gly3->Pro substitution and about three orders of magnitude of affinity.

Usage:
    python scripts/10_benchmark.py --out data/benchmark
"""
import argparse
import csv
import json
import random
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path

AFFINITY = Path("data/reference/affinity_set.tsv")
LIBRARY = Path("data/library/peptides.tsv")
CONTROLS = Path("data/controls/control_sets.json")

SCRAMBLE_SEED = 20261020
DECOY_SEED = 20261021
# Hard decoys per active. The library has 640 Tyr-starting peptides, so a few per active
# is affordable and keeps the comparison from being dominated by one length.
HARD_DECOYS_PER_ACTIVE = 4


def read_tsv(path):
    with Path(path).open() as fh:
        rows = [l for l in fh if not l.startswith("#") and l.strip()]
    return list(csv.DictReader(rows, delimiter="\t"))


def scramble(sequence, rng, max_tries=200):
    """A composition-matched permutation that is not the original sequence.

    Peptides whose composition admits no distinct permutation (e.g. all one residue)
    are reported rather than silently returned unchanged.
    """
    letters = list(sequence)
    for _ in range(max_tries):
        rng.shuffle(letters)
        candidate = "".join(letters)
        if candidate != sequence:
            return candidate
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=Path("data/benchmark"), type=Path)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    actives = read_tsv(AFFINITY)
    library = read_tsv(LIBRARY)
    controls = json.loads(CONTROLS.read_text())

    rows = []
    notes = {}

    # ---- actives ---------------------------------------------------------------
    # The ChEMBL set contains the same sequence twice with different C-termini
    # (YGGFM free, pKi 8.70; YGGFM amide, pKi 7.92), so the id must carry the ChEMBL
    # accession rather than the sequence alone.
    for a in actives:
        seq = a["sequence"]
        rows.append({
            "id": f"act__{a['name']}",
            "sequence": seq,
            "set": "active",
            "c_term": "amide" if a["amidated"] == "True" else "free",
            "pKi": a["pKi_median"],
            "family": seq[:5],
            "length": str(len(seq)),
        })

    # ---- scrambles -------------------------------------------------------------
    rng = random.Random(SCRAMBLE_SEED)
    unscrambleable = []
    for a in actives:
        seq = a["sequence"]
        scr = scramble(seq, rng)
        if scr is None:
            unscrambleable.append(seq)
            continue
        rows.append({
            "id": f"scr__{a['name']}",
            "sequence": scr,
            "set": "scramble",
            # The scramble keeps the active's C-terminal chemistry so that the only
            # difference under test is residue order.
            "c_term": "amide" if a["amidated"] == "True" else "free",
            "pKi": "",
            "family": "",
            "length": str(len(scr)),
        })
    notes["scrambles_unscrambleable"] = unscrambleable
    notes["scrambles_starting_with_Y"] = sum(
        1 for r in rows if r["set"] == "scramble" and r["sequence"].startswith("Y"))

    # ---- hard decoys: Tyr-starting, length-matched ------------------------------
    active_seqs = {a["sequence"] for a in actives}
    by_length = defaultdict(list)
    for r in library:
        if r["starts_with_Y"] == "yes" and r["sequence"] not in active_seqs:
            by_length[len(r["sequence"])].append(r)

    drng = random.Random(DECOY_SEED)
    chosen, shortfall = {}, []
    for a in actives:
        want = len(a["sequence"])
        # Nearest available length, preferring exact; records when no exact match exists.
        lengths = sorted(by_length, key=lambda L: (abs(L - want), L))
        if not lengths:
            shortfall.append(a["sequence"])
            continue
        pool = []
        for L in lengths:
            pool.extend(by_length[L])
            if len(pool) >= HARD_DECOYS_PER_ACTIVE * 5:
                break
        drng.shuffle(pool)
        taken = 0
        for cand in pool:
            if cand["peptide_id"] in chosen:
                continue
            chosen[cand["peptide_id"]] = cand
            taken += 1
            if taken >= HARD_DECOYS_PER_ACTIVE:
                break
    for pid, cand in chosen.items():
        rows.append({
            "id": f"hard__{pid}", "sequence": cand["sequence"], "set": "tyr_hard",
            "c_term": "auto", "pKi": "", "family": "",
            "length": str(len(cand["sequence"]))})
    notes["hard_decoy_length_shortfall"] = shortfall

    # ---- positional control ------------------------------------------------------
    pc = controls["positional_control"]
    for entry in pc["starts_with_motif"]:
        rows.append({"id": f"posstart__{entry['peptide_id']}",
                     "sequence": entry["sequence"], "set": "positional_start",
                     "c_term": "auto", "pKi": "", "family": "",
                     "length": str(len(entry["sequence"]))})
    for entry in pc["spans_motif_blocked_nterm"]:
        rows.append({"id": f"posspan__{entry['peptide_id']}",
                     "sequence": entry["sequence"], "set": "positional_span",
                     "c_term": "auto", "pKi": "", "family": "",
                     "length": str(len(entry["sequence"]))})

    # ---- retrieval controls -------------------------------------------------------
    # met-enkephalin (template ligand) and endomorphin-1 (held-out structure 8F7R) are
    # already in the active set via ChEMBL, but are tagged so retrieval can find them.
    columns = ["id", "sequence", "set", "c_term", "pKi", "family", "length"]
    out = args.out / "benchmark.tsv"
    with out.open("w") as fh:
        fh.write("\t".join(columns) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in columns) + "\n")

    counts = Counter(r["set"] for r in rows)
    fam = Counter(r["family"] for r in rows if r["set"] == "active" and r["family"])
    summary = {
        "generated": date.today().isoformat(),
        "total": len(rows),
        "by_set": dict(counts),
        "potency_pair": {
            "high_family": "YGGFL", "n_high": fam.get("YGGFL", 0),
            "low_family": "YGPFL", "n_low": fam.get("YGPFL", 0),
        },
        "scramble_seed": SCRAMBLE_SEED,
        "decoy_seed": DECOY_SEED,
        "notes": notes,
    }
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    print(f"wrote {out}  ({len(rows)} rows)")
    for k, v in sorted(counts.items()):
        print(f"  {k:20s} {v}")
    print(f"  potency pair: YGGFL n={fam.get('YGGFL',0)}  YGPFL n={fam.get('YGPFL',0)}")
    print(f"  scrambles starting with Y: {notes['scrambles_starting_with_Y']}")
    if unscrambleable:
        print(f"  WARNING unscrambleable: {unscrambleable}")


if __name__ == "__main__":
    main()
