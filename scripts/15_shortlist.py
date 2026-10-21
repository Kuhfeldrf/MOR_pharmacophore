#!/usr/bin/env python3
"""The library shortlist, ranked within the strata the actives actually occupy.

Post-hoc, like 14. The reason this script exists instead of the raw tiers from 13:

The score is close to a staircase in the number of aromatic residues (F/W/Y/H) in the
scored window - adjacent strata barely overlap, and Spearman rho between score and
aromatic count over all 14,571 library windows is +0.909. That is partly real signal:
the template's two anchor features are the Tyr1 and Phe4 aromatics of the opioid
message, and a window with no aromatic group cannot present them. But it also means a
raw cross-library ranking is mostly a count, and the top of it is a His-rich repeat
family in strata of 19 and 5 windows that contain no known active.

Full percentile stratification was tried (script 14) and REJECTED: it failed its own
control (E1 moved 0.803 -> 0.712 when the aromatic-matched family contrast should have
been invariant, because family members fall in different strata at window 8), and it
promoted zero-aromatic peptides to the 99.9th percentile for being the best of a
stratum in which nothing matches the template at all.

What is reported instead is narrower and survives both objections: rank within the
aromatic strata that the endogenous ligands occupy (2 and 3), against a bar set by the
weakest endogenous opioid peptide in the same stratum. Peptides in the 4- and 5-aromatic
strata are listed separately and flagged, not ranked alongside - 24 windows is not a
calibratable background.

Usage:
    python scripts/15_shortlist.py --library results/library_w8 --out docs
"""
import argparse
import csv
import glob
import json
from collections import defaultdict
from datetime import date
from pathlib import Path

SCORE = "tversky_ref_mean"
AROMATIC = "FWYH"
# Strata containing at least one endogenous ligand; others cannot be calibrated here.
RANKABLE_STRATA = (2, 3)
MIN_STRATUM_WINDOWS = 100


def n_aromatic(window):
    return sum(window.count(c) for c in AROMATIC)


def load_blocks(d):
    out = []
    for f in sorted(glob.glob(str(Path(d) / "block*.tsv"))):
        with open(f) as fh:
            out.extend(csv.DictReader(fh, delimiter="\t"))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--library", default=Path("results/library_w8"), type=Path)
    ap.add_argument("--out", default=Path("docs"), type=Path)
    ap.add_argument("--top", default=40, type=int)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    rows = load_blocks(args.library)
    lib = [r for r in rows if r["set"] == "library"]
    refs = [r for r in rows if r["set"] != "library"]

    # ---- the bar, per stratum: the weakest ENDOGENOUS ligand in that stratum --------
    endo = defaultdict(list)
    for r in refs:
        if r["set"] in ("positive_endogenous", "template_ligand"):
            endo[n_aromatic(r["window"])].append(
                (float(r[SCORE]), r["id"].replace("ref__", "")))
    bars = {k: min(v) for k, v in endo.items()}

    # ---- stratum backgrounds (distinct windows) ------------------------------------
    bg = defaultdict(list)
    seen = set()
    for r in lib:
        if r["window"] in seen:
            continue
        seen.add(r["window"])
        bg[n_aromatic(r["window"])].append(float(r[SCORE]))
    for k in bg:
        bg[k].sort()

    shortlist, flagged = [], []
    for r in lib:
        a = n_aromatic(r["window"])
        v = float(r[SCORE])
        pop = bg[a]
        rec = {
            "peptide_id": r["id"], "sequence": r["sequence"], "window": r["window"],
            "n_aromatic": a, "score": v,
            "stratum_percentile": sum(1 for x in pop if x <= v) / len(pop),
            "free_nterm_tyr": r["sequence"].startswith("Y"),
            "truncated": r["truncated"] == "True",
            "matched_features": int(float(r["matched_features"])),
        }
        if a in RANKABLE_STRATA and a in bars and v >= bars[a][0]:
            rec["bar"] = bars[a][1]
            rec["bar_score"] = bars[a][0]
            shortlist.append(rec)
        elif a > max(RANKABLE_STRATA):
            rec["flag"] = ("stratum of %d windows: too thin to calibrate, no known "
                           "active in it" % len(pop))
            flagged.append(rec)
    shortlist.sort(key=lambda r: (-r["n_aromatic"], -r["score"]))
    flagged.sort(key=lambda r: -r["score"])

    def tsv(path, recs, extra):
        cols = (["rank", "peptide_id", "sequence", "window", "n_aromatic", "score",
                 "stratum_percentile", "free_nterm_tyr", "matched_features"] + extra)
        with path.open("w") as fh:
            fh.write("\t".join(cols) + "\n")
            for i, r in enumerate(recs, 1):
                r = dict(r, rank=i)
                fh.write("\t".join(
                    ("%.4f" % r[c]) if c in ("score", "stratum_percentile", "bar_score")
                    else str(r.get(c, "")) for c in cols) + "\n")

    tsv(args.out / "library_shortlist.tsv", shortlist, ["bar", "bar_score"])
    tsv(args.out / "library_high_aromatic_flagged.tsv", flagged, ["flag"])

    by_stratum = {}
    for a in sorted(bg):
        v = bg[a]
        by_stratum[a] = {
            "distinct_windows": len(v), "median": v[len(v) // 2], "max": v[-1],
            "calibratable": len(v) >= MIN_STRATUM_WINDOWS,
            "endogenous_bar": bars.get(a, [None, None])[0],
            "endogenous_bar_ligand": bars.get(a, [None, None])[1],
            "n_shortlisted": sum(1 for r in shortlist if r["n_aromatic"] == a),
        }

    payload = {
        "generated": date.today().isoformat(), "post_hoc": True,
        "score_field": SCORE,
        "rule": ("rank within aromatic strata occupied by endogenous ligands "
                 "(%s); bar = weakest endogenous ligand in the same stratum"
                 % list(RANKABLE_STRATA)),
        "rejected_alternative": ("full percentile stratification (script 14) - failed "
                                 "its own E1 invariance control and promoted "
                                 "non-aromatic peptides"),
        "strata": by_stratum,
        "n_shortlisted": len(shortlist), "n_flagged_high_aromatic": len(flagged),
        "n_shortlisted_tyr_start": sum(1 for r in shortlist if r["free_nterm_tyr"]),
        "shortlist_top": shortlist[:args.top],
        "flagged_top": flagged[:15],
    }
    (args.out / "library_shortlist.json").write_text(json.dumps(payload, indent=2) + "\n")

    print("Aromatic strata of the library (distinct windows):")
    print("  arom  windows   median     max    bar (weakest endogenous)   shortlisted")
    for a, d in by_stratum.items():
        bar = ("%.4f %s" % (d["endogenous_bar"], d["endogenous_bar_ligand"])
               if d["endogenous_bar"] else "-")
        print("  %4d  %7d   %.4f  %.4f   %-26s %s" % (
            a, d["distinct_windows"], d["median"], d["max"], bar,
            d["n_shortlisted"] or "-"))
    print("\nshortlist: %d peptides (%d Tyr-starting) in strata %s"
          % (len(shortlist), payload["n_shortlisted_tyr_start"],
             list(RANKABLE_STRATA)))
    print("flagged high-aromatic: %d peptides, reported separately" % len(flagged))
    print("\nTop %d shortlist:" % min(args.top, len(shortlist)))
    print("  rank sequence                        win       arom  score   pct   Y")
    for i, r in enumerate(shortlist[:args.top], 1):
        print("  %4d %-30s %-9s %d    %.4f %.3f  %s" % (
            i, r["sequence"][:30], r["window"], r["n_aromatic"], r["score"],
            r["stratum_percentile"], "Y" if r["free_nterm_tyr"] else ""))
    print("\nwrote %s and %s" % (args.out / "library_shortlist.tsv",
                                 args.out / "library_high_aromatic_flagged.tsv"))


if __name__ == "__main__":
    main()
