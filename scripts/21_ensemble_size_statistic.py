#!/usr/bin/env python3
"""Is the published ensemble-size penalty specific to the `max` statistic? (D4.2)

McNutt et al. (J Chem Inf Model 2023, 63:6598) recommend fewer than 25 conformers for
pharmacophore screening, and are explicit that the reason is accuracy rather than cost:
larger ensembles "increase the false positive rate faster than the true positive rate".
This pipeline uses 800. Taken at face value, its most expensive decision is wrong by a
factor of thirty.

The proposed reconciliation (DECISIONS.md D4.2) is mechanistic. Conventional
pharmacophore screening is hit-or-miss - one conformer matching is enough - which is the
`max` statistic. Under `max`, every extra conformer is another chance at a spurious
match, so false positives accumulate with ensemble size. Under `mean`, extra conformers
refine an average instead of adding draws to a maximum, so the same mechanism does not
apply.

That predicts an interaction, which is what this script measures: going from 400 to 800
conformers at window 8, discrimination should DEGRADE (or at best not improve) under
`max`, and hold or improve under `mean` and `top10`.

If instead `mean` degrades too, the reconciliation fails and the ensemble size is simply
too large, D4.2 is wrong, and the production configuration needs revisiting.

Pure re-analysis: both conformer counts, both seeds, all statistics are recomputable from
score tables already on disk.

Usage:
    python scripts/21_ensemble_size_statistic.py --stage results/stageB --out docs
"""
import argparse
import csv
import glob
import json
import re
from collections import defaultdict
from datetime import date
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from mor_pharmacophore.alignit import read_scores  # noqa: E402

METRIC = "tversky_ref"
STATISTICS = ("max", "mean", "top10")
AROMATIC = "FWYH"


def auroc(pos, neg):
    if not pos or not neg:
        return float("nan")
    wins = sum(1.0 if p > n else (0.5 if p == n else 0.0) for p in pos for n in neg)
    return wins / (len(pos) * len(neg))


def _ranks(xs):
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    r = [0.0] * len(xs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        for k in range(i, j + 1):
            r[order[k]] = (i + j) / 2.0 + 1.0
        i = j + 1
    return r


def spearman(xs, ys):
    if len(xs) < 3:
        return float("nan")
    rx, ry = _ranks(xs), _ranks(ys)
    n = len(xs)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = sum((a - mx) ** 2 for a in rx) ** 0.5
    dy = sum((b - my) ** 2 for b in ry) ** 0.5
    return num / (dx * dy) if dx and dy else float("nan")


def collapse(rows, statistic):
    g = defaultdict(list)
    for r in rows:
        g[r["candidate"]].append(r[METRIC])
    out = {}
    for name, vals in g.items():
        vals.sort(reverse=True)
        if statistic == "max":
            out[name] = vals[0]
        elif statistic == "mean":
            out[name] = sum(vals) / len(vals)
        else:
            k = max(1, len(vals) // 10)
            out[name] = sum(vals[:k]) / k
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default=Path("results/stageB"), type=Path)
    ap.add_argument("--benchmark", default=Path("data/benchmark/benchmark.tsv"),
                    type=Path)
    ap.add_argument("--window", default=8, type=int)
    ap.add_argument("--out", default=Path("docs"), type=Path)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    with args.benchmark.open() as fh:
        bench = {r["id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    sets, fam = defaultdict(list), defaultdict(list)
    for bid, r in bench.items():
        sets[r["set"]].append(bid)
        if r["set"] == "active" and r["family"]:
            fam[r["family"]].append(bid)

    runs = defaultdict(dict)
    for d in sorted(args.stage.glob("nc*_w*_s*")):
        m = re.match(r"nc(\d+)_w(\d+)_s(\d+)", d.name)
        if not m:
            continue
        nc, w, seed = int(m.group(1)), int(m.group(2)), int(m.group(3))
        if w != args.window:
            continue
        rows = []
        for tab in sorted(glob.glob(str(d / "block*" / "scores.tab"))):
            rows.extend(read_scores(tab))
        if rows:
            runs[nc][seed] = rows

    results = []
    for nc in sorted(runs):
        seeds = sorted(runs[nc])
        if len(seeds) < 2:
            continue
        for stat in STATISTICS:
            a = collapse(runs[nc][seeds[0]], stat)
            b = collapse(runs[nc][seeds[1]], stat)
            shared = sorted(set(a) & set(b))

            def g(ids, t):
                return [t[i] for i in ids if i in t]

            # averaged over the two seeds, so a single lucky seed cannot drive the trend
            s2 = sum(auroc(g(sets["active"], t), g(sets["scramble"], t))
                     for t in (a, b)) / 2
            e3 = sum(auroc(g(sets["active"], t), g(sets["tyr_hard"], t))
                     for t in (a, b)) / 2
            e1 = sum(auroc(g(fam.get("YGGFL", []), t), g(fam.get("YGPFL", []), t))
                     for t in (a, b)) / 2
            results.append({
                "n_conformers": nc, "statistic": stat,
                "S1_spearman": spearman([a[k] for k in shared],
                                        [b[k] for k in shared]),
                "S2_auroc_vs_scramble": s2,
                "E3_auroc_vs_tyr_hard": e3,
                "E1_potency_pair_auroc": e1,
                "n_scored": len(shared),
            })

    # ---- the interaction: change from the smallest to the largest ensemble --------
    counts = sorted({r["n_conformers"] for r in results})
    deltas = {}
    if len(counts) >= 2:
        lo, hi = counts[0], counts[-1]
        for stat in STATISTICS:
            rl = next((r for r in results
                       if r["n_conformers"] == lo and r["statistic"] == stat), None)
            rh = next((r for r in results
                       if r["n_conformers"] == hi and r["statistic"] == stat), None)
            if rl and rh:
                deltas[stat] = {
                    "from": lo, "to": hi,
                    "dS2": rh["S2_auroc_vs_scramble"] - rl["S2_auroc_vs_scramble"],
                    "dE3": rh["E3_auroc_vs_tyr_hard"] - rl["E3_auroc_vs_tyr_hard"],
                    "dS1": rh["S1_spearman"] - rl["S1_spearman"],
                }

    verdict = None
    if "max" in deltas and "mean" in deltas:
        max_degrades = deltas["max"]["dS2"] < 0 or deltas["max"]["dE3"] < 0
        mean_holds = deltas["mean"]["dS2"] >= -0.01 and deltas["mean"]["dE3"] >= -0.01
        verdict = ("reconciliation supported: the penalty is max-specific"
                   if (max_degrades and mean_holds)
                   else "reconciliation NOT supported")

    payload = {
        "generated": date.today().isoformat(), "metric": METRIC,
        "window": args.window,
        "hypothesis": ("McNutt 2023's larger-ensemble penalty is a false-positive "
                       "mechanism specific to hit-or-miss (max) matching, and should "
                       "not apply to the ensemble mean"),
        "results": results, "deltas_small_to_large_ensemble": deltas,
        "verdict": verdict,
        "caveat": ("Two conformer counts and one window is a weak test of an "
                   "interaction; it can refute the reconciliation but cannot establish "
                   "it. AUROC is also the wrong summary for screening (D6.8) - used "
                   "here only because the question is a relative change in the same "
                   "metric across ensemble sizes."),
    }
    (args.out / "ensemble_size_statistic.json").write_text(
        json.dumps(payload, indent=2) + "\n")

    print("Ensemble size x ensemble statistic, window %d, %s" % (args.window, METRIC))
    print("(S2 and E3 averaged over both seeds)\n")
    print("  nconf  stat     S1rho    S2      E3      E1")
    for r in sorted(results, key=lambda r: (r["statistic"], r["n_conformers"])):
        print("  %5d  %-6s  %6.3f  %6.3f  %6.3f  %6.3f" % (
            r["n_conformers"], r["statistic"], r["S1_spearman"],
            r["S2_auroc_vs_scramble"], r["E3_auroc_vs_tyr_hard"],
            r["E1_potency_pair_auroc"]))
    if deltas:
        print("\n  change from %d to %d conformers:" % (counts[0], counts[-1]))
        print("  stat     dS1      dS2      dE3")
        for stat, d in deltas.items():
            print("  %-6s  %+6.3f  %+6.3f  %+6.3f"
                  % (stat, d["dS1"], d["dS2"], d["dE3"]))
    if verdict:
        print("\n  %s" % verdict)
    print("\nwrote %s" % (args.out / "ensemble_size_statistic.json"))


if __name__ == "__main__":
    main()
