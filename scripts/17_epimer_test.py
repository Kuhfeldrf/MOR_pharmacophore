#!/usr/bin/env python3
"""E7: epimer discrimination - does the score read shape, or count aromatics?

This is the test the study should have been built around. Every other contrast available
here varies composition along with shape, and docs/optimisation_study.md 7 showed that
across the library the score tracks aromatic-residue count at rho +0.909, while
scripts/16_reselect_with_e6.py showed that no choice of metric separates the two: every
metric that weakens the confound weakens E1 in proportion, because covering the
template's aromatic anchors IS what the score measures.

An L/D epimer pair breaks that. The two molecules have identical atoms, identical
aromatic count, identical formula and identical constitution; they differ only in the
handedness of one alpha carbon. No composition rule of any kind can tell them apart, so
any separation the screen achieves is necessarily conformational. If the method cannot
order epimer pairs, then what E1 measured was the composition difference that happens to
accompany Gly3->Pro, and the pharmacophore matching is adding nothing over a feature
count.

Direction and what counts as a pass
-----------------------------------
D-Ala2 is required for mu affinity in the dermorphin family; the L-Ala2 diastereomer is
markedly weaker. The literature evidence for the direction is qualitative rather than a
matched pKi pair, so this is scored as a SIGN TEST over pairs - does score(D) > score(L)
- not as an AUROC against measured affinity. With n pairs the null is a coin flip, so
the exact binomial p-value is reported and no confidence interval is invented for a
quantity that has no measured scale here.

Pair selection caveat, stated before the numbers: the dermorphin and [D-Ala2]enkephalin
amide pairs are mu-relevant and carry the directional prediction. The deltorphins are
DELTA-selective - their D-Ala2 confers delta preference - so against a mu template there
is no clear prediction for them, and they are reported separately and excluded from the
headline sign test rather than counted in whichever direction happens to help.

Usage:
    python scripts/17_epimer_test.py --run results/epimers \\
        --pairs data/epimers/epimer_pairs.tsv --out docs
"""
import argparse
import csv
import glob
import json
from collections import defaultdict
from datetime import date
from math import comb
from pathlib import Path

METRIC = "tversky_ref_mean"
# mu-relevant families carrying the directional prediction
MU_FAMILIES = ("dermorphin", "dala2_met_enkephalin", "dala2_leu_enkephalin")
AROMATIC = "FWYH"


def binom_p_two_sided(k, n):
    """Exact two-sided binomial p-value against p=0.5."""
    if n == 0:
        return float("nan")
    probs = [comb(n, i) / 2 ** n for i in range(n + 1)]
    obs = probs[k]
    return min(1.0, sum(p for p in probs if p <= obs + 1e-12))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default=Path("results/epimers"), type=Path)
    ap.add_argument("--pairs", default=Path("data/epimers/epimer_pairs.tsv"), type=Path)
    ap.add_argument("--out", default=Path("docs"), type=Path)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    with args.pairs.open() as fh:
        meta = {r["id"]: r for r in csv.DictReader(fh, delimiter="\t")}

    scored = {}
    for f in sorted(glob.glob(str(args.run / "block*.tsv"))) or \
            sorted(glob.glob(str(args.run / "*.tsv"))):
        with open(f) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                scored[r["id"]] = r
    if not scored:
        raise SystemExit(f"no score rows found under {args.run}")

    pairs = defaultdict(dict)
    for pid, m in meta.items():
        if pid in scored:
            pairs[m["pair"]][m["form"]] = (scored[pid], m)

    rows = []
    for pair, forms in sorted(pairs.items()):
        if "D" not in forms or "L" not in forms:
            continue
        d_row, d_meta = forms["D"]
        l_row, l_meta = forms["L"]
        ds, ls = float(d_row[METRIC]), float(l_row[METRIC])
        ard = sum(d_row["window"].count(c) for c in AROMATIC)
        arl = sum(l_row["window"].count(c) for c in AROMATIC)
        rows.append({
            "pair": pair,
            "mu_relevant": pair in MU_FAMILIES,
            "D_sequence": d_meta["sequence"], "L_sequence": l_meta["sequence"],
            "D_score": ds, "L_score": ls, "delta": ds - ls,
            "correct_direction": ds > ls,
            "aromatic_count_D": ard, "aromatic_count_L": arl,
            "aromatic_matched": ard == arl,
            "D_conformers": int(d_row["n_conformers"]),
            "L_conformers": int(l_row["n_conformers"]),
            "D_matched_features": int(float(d_row["matched_features"])),
            "L_matched_features": int(float(l_row["matched_features"])),
        })

    mu = [r for r in rows if r["mu_relevant"]]
    k = sum(1 for r in mu if r["correct_direction"])
    n = len(mu)
    delta_tracked = [r for r in rows if not r["aromatic_matched"]]

    payload = {
        "generated": date.today().isoformat(),
        "metric": METRIC,
        "configuration": {"window": 8, "n_conformers": 800},
        "test": "sign test over mu-relevant epimer pairs, exact two-sided binomial",
        "n_mu_pairs": n, "n_correct": k,
        "binomial_p_two_sided": binom_p_two_sided(k, n),
        "all_pairs_aromatic_matched": all(r["aromatic_matched"] for r in rows),
        "pairs": rows,
        "note": ("Delta-selective deltorphin pairs are reported but excluded from the "
                 "sign test: against a mu template they carry no directional "
                 "prediction."),
    }
    (args.out / "epimer_test.json").write_text(json.dumps(payload, indent=2) + "\n")

    print("E7 epimer discrimination, window 8 / 800 conformers, %s\n" % METRIC)
    print("  pair                      D-seq      L-seq      D score  L score   delta"
          "   arom  dir")
    for r in sorted(rows, key=lambda r: (not r["mu_relevant"], -r["delta"])):
        print("  %-24s %-10s %-10s %.4f   %.4f  %+.4f  %d=%d  %s%s" % (
            r["pair"], r["D_sequence"], r["L_sequence"], r["D_score"], r["L_score"],
            r["delta"], r["aromatic_count_D"], r["aromatic_count_L"],
            "OK " if r["correct_direction"] else "NO ",
            "" if r["mu_relevant"] else " (delta-selective, excluded)"))
    print("\n  every pair aromatic-matched: %s"
          % payload["all_pairs_aromatic_matched"])
    if delta_tracked:
        print("  WARNING pairs not aromatic-matched: %s"
              % [r["pair"] for r in delta_tracked])
    print("\n  mu-relevant sign test: %d/%d in the predicted direction, "
          "exact two-sided p = %.4f" % (k, n, binom_p_two_sided(k, n)))
    if n and k == n:
        print("  All mu pairs ordered correctly at exactly matched composition.")
    elif n and k <= n / 2:
        print("  The screen does not order epimer pairs: at exactly matched feature\n"
              "  composition it has no resolving power. Note what this does and does "
              "not\n  say about E1. Gly3->Pro is aromatic-matched, so E1 is not "
              "explained by the\n  aromatic confound - but Pro adds a lipophilic ring "
              "and removes a backbone\n  donor, so it still changes FEATURE "
              "composition. Epimerisation changes no\n  feature count at all. The "
              "pattern across the three is that the score tracks\n  pharmacophore "
              "feature composition and not 3D arrangement; see\n  "
              "scripts/18_effect_sizes.py for each contrast's margin in jitter units.")
    print("\nwrote %s" % (args.out / "epimer_test.json"))


if __name__ == "__main__":
    main()
