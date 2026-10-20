#!/usr/bin/env python3
"""Stage A: choose the configuration using the SELECTION metrics only (STUDY.md 3).

This script deliberately computes nothing else. It is given the selection set
(actives + composition-matched scrambles) and can only compute:

  S1 reproducibility  Spearman rho between two independent conformer seeds (gate >= 0.95)
  S2 classification   AUROC, actives vs scrambles

The evaluation metrics - the potency pair, potency ranking, added value over the Tyr
baseline, held-out retrieval and the positional rule - need peptides that are not in
this set, so they cannot be computed here even by accident. That is the point: the
separation between selecting and evaluating is enforced by what the input contains,
not only by intent.

Usage:
    python scripts/11_select.py --stage results/stageA --out docs
"""
import argparse
import csv
import json
import re
from collections import defaultdict
from datetime import date
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from mor_pharmacophore.alignit import read_scores  # noqa: E402

METRICS = ("tanimoto", "tversky_ref", "tversky_db")
STATISTICS = ("max", "mean", "top10")
S1_GATE = 0.95


def auroc(pos, neg):
    if not pos or not neg:
        return float("nan")
    wins = sum(1.0 if p > n else (0.5 if p == n else 0.0) for p in pos for n in neg)
    return wins / (len(pos) * len(neg))


def _ranks(xs):
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    ranks = [0.0] * len(xs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        for k in range(i, j + 1):
            ranks[order[k]] = (i + j) / 2.0 + 1.0
        i = j + 1
    return ranks


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


def collapse(rows, metric, statistic):
    grouped = defaultdict(list)
    for r in rows:
        grouped[r["candidate"]].append(r[metric])
    out = {}
    for name, vals in grouped.items():
        vals.sort(reverse=True)
        if statistic == "max":
            out[name] = vals[0]
        elif statistic == "mean":
            out[name] = sum(vals) / len(vals)
        else:
            k = max(1, len(vals) // 10)
            out[name] = sum(vals[:k]) / k
    return out


def load_run(run_dir):
    """Concatenate every block's per-conformer score table for one configuration."""
    rows = []
    for tab in sorted(Path(run_dir).glob("block*/scores.tab")):
        rows.extend(read_scores(tab))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default=Path("results/stageA"), type=Path)
    ap.add_argument("--selection", default=Path("data/benchmark/selection_set.tsv"),
                    type=Path)
    ap.add_argument("--out", default=Path("docs"), type=Path)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    with args.selection.open() as fh:
        bench = {r["id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    actives = [i for i, r in bench.items() if r["set"] == "active"]
    scrambles = [i for i, r in bench.items() if r["set"] == "scramble"]
    assert actives and scrambles, "selection set must contain actives and scrambles"

    runs = defaultdict(dict)
    for d in sorted(args.stage.glob("nc*_w*_s*")):
        m = re.match(r"nc(\d+)_w(\d+)_s(\d+)", d.name)
        if not m:
            continue
        nc, w, seed = int(m.group(1)), int(m.group(2)), int(m.group(3))
        rows = load_run(d)
        if rows:
            runs[(nc, w)][seed] = rows

    results = []
    for (nc, w), by_seed in sorted(runs.items()):
        seeds = sorted(by_seed)
        if len(seeds) < 2:
            print(f"  skipping nc={nc} w={w}: only {len(seeds)} seed(s)")
            continue
        for metric in METRICS:
            for stat in STATISTICS:
                a = collapse(by_seed[seeds[0]], metric, stat)
                b = collapse(by_seed[seeds[1]], metric, stat)
                shared = sorted(set(a) & set(b))
                if len(shared) < 10:
                    continue
                rho = spearman([a[k] for k in shared], [b[k] for k in shared])
                deltas = sorted(abs(a[k] - b[k]) for k in shared)
                vals = sorted(a[k] for k in shared)
                spread = vals[int(0.9 * len(vals))] - vals[int(0.1 * len(vals))]
                s2 = auroc([a[i] for i in actives if i in a],
                           [a[i] for i in scrambles if i in a])
                s2b = auroc([b[i] for i in actives if i in b],
                            [b[i] for i in scrambles if i in b])
                results.append({
                    "n_conformers": nc, "window": w, "metric": metric,
                    "statistic": stat, "n_scored": len(shared),
                    "S1_spearman": rho,
                    "S1_median_abs_delta": deltas[len(deltas) // 2],
                    "S1_p10_p90_spread": spread,
                    "S1_delta_over_spread": (deltas[len(deltas) // 2] / spread
                                             if spread else float("nan")),
                    "S1_pass": bool(rho == rho and rho >= S1_GATE),
                    "S2_auroc_seed1": s2, "S2_auroc_seed2": s2b,
                    "S2_auroc_mean": (s2 + s2b) / 2,
                })

    eligible = [r for r in results if r["S1_pass"]]
    chosen = max(eligible, key=lambda r: (r["S2_auroc_mean"], -r["n_conformers"])) \
        if eligible else None

    # STUDY.md 4.1: the window decides which questions are answerable, so the selection
    # metrics choose n_conformers/metric/statistic WITHIN each window and all windows are
    # carried into evaluation.
    per_window = {}
    for w in sorted({r["window"] for r in results}):
        pool = [r for r in eligible if r["window"] == w]
        if pool:
            per_window[str(w)] = max(
                pool, key=lambda r: (r["S2_auroc_mean"], -r["n_conformers"]))

    payload = {"generated": date.today().isoformat(), "s1_gate": S1_GATE,
               "selection_metrics_only": True,
               "n_configurations": len(results),
               "n_passing_S1": len(eligible),
               "chosen": chosen, "chosen_per_window": per_window,
               "results": results}
    (args.out / "selection_results.json").write_text(json.dumps(payload, indent=2) + "\n")

    print(f"{len(results)} configurations; {len(eligible)} pass S1 (rho >= {S1_GATE})\n")
    print("  nconf win metric      stat    S1rho  d/spread  S2auc  gate")
    for r in sorted(results, key=lambda r: (-r["S1_spearman"], -r["S2_auroc_mean"]))[:18]:
        print("  %5d %3d %-11s %-6s %6.3f %8.3f %6.3f  %s" % (
            r["n_conformers"], r["window"], r["metric"], r["statistic"],
            r["S1_spearman"], r["S1_delta_over_spread"], r["S2_auroc_mean"],
            "PASS" if r["S1_pass"] else ""))
    print("\nBest per window (S1 gate + S2 only):")
    for w, r in sorted(per_window.items(), key=lambda kv: int(kv[0])):
        print("  window %-3s nconf=%-4d %-11s %-6s  S1rho=%.3f  S2=%.3f" % (
            w, r["n_conformers"], r["metric"], r["statistic"],
            r["S1_spearman"], r["S2_auroc_mean"]))
    if chosen:
        print("\nCHOSEN on S1 gate + S2 only:")
        for k in ("n_conformers", "window", "metric", "statistic",
                  "S1_spearman", "S2_auroc_mean"):
            print(f"  {k:18s} {chosen[k]}")
    else:
        print("\nNO configuration passed the S1 reproducibility gate.")


if __name__ == "__main__":
    main()
