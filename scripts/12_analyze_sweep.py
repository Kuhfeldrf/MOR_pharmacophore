#!/usr/bin/env python3
"""Analyse the optimisation sweep (STUDY.md 3-4).

Selection metrics (may choose the configuration):
  S1 reproducibility  Spearman rho between two conformer seeds; gate rho >= 0.95
  S2 classification   AUROC, actives vs composition-matched scrambles

Evaluation metrics (reported, never used to choose):
  E1 potency pair     AUROC, YGGFL family (pKi 7.4-8.7) vs YGPFL family (4.7-5.5)
  E2 potency ranking  Spearman rho vs pKi over the actives
  E3 added value      AUROC, actives vs Tyr-starting length-matched library decoys
  E4 held-out ligand  endomorphin-1 percentile against the decoy background
  E5 positional rule  AUROC, starts-with-YPFP vs spans-YPFP

The split is enforced here in code, not only in prose: `choose_configuration` is given
the selection metrics alone.

Usage:
    python scripts/12_analyze_sweep.py --sweep results/sweep --out docs
"""
import argparse
import csv
import json
import random
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


# ---------------------------------------------------------------- statistics
def auroc(pos, neg):
    if not pos or not neg:
        return float("nan")
    wins = 0.0
    for p in pos:
        for n in neg:
            wins += 1.0 if p > n else (0.5 if p == n else 0.0)
    return wins / (len(pos) * len(neg))


def bootstrap_auroc(pos, neg, n=2000, seed=7):
    if not pos or not neg:
        return float("nan"), float("nan")
    rng = random.Random(seed)
    vals = sorted(auroc([rng.choice(pos) for _ in pos],
                        [rng.choice(neg) for _ in neg]) for _ in range(n))
    return vals[int(0.025 * n)], vals[int(0.975 * n) - 1]


def _ranks(xs):
    """Average ranks, so ties are handled correctly."""
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    ranks = [0.0] * len(xs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = avg
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


def percentile_of(value, population):
    if not population:
        return float("nan")
    return sum(1 for p in population if p <= value) / len(population)


# ---------------------------------------------------------------- loading
def collapse(rows, metric, statistic):
    """One score per candidate, under the chosen metric and ensemble statistic."""
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
        else:                                  # top10: mean of the best decile
            k = max(1, len(vals) // 10)
            out[name] = sum(vals[:k]) / k
    return out


def load_config(path):
    rows = []
    for tab in sorted(Path(path).glob("block*/scores.tab")):
        rows.extend(read_scores(tab))
    return rows or None


def parse_tag(tag):
    m = re.match(r"nc(\d+)_w(\d+)_s(\d+)", tag)
    return (int(m.group(1)), int(m.group(2)), int(m.group(3))) if m else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sweep", default=Path("results/sweep"), type=Path)
    ap.add_argument("--benchmark", default=Path("data/benchmark/benchmark.tsv"),
                    type=Path)
    ap.add_argument("--out", default=Path("docs"), type=Path)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    with args.benchmark.open() as fh:
        bench = {r["id"]: r for r in csv.DictReader(fh, delimiter="\t")}

    sets = defaultdict(list)
    for bid, r in bench.items():
        sets[r["set"]].append(bid)
    fam = defaultdict(list)
    for bid, r in bench.items():
        if r["set"] == "active" and r["family"]:
            fam[r["family"]].append(bid)

    # group runs by (n_conformers, window), keeping both seeds
    runs = defaultdict(dict)
    for d in sorted(args.sweep.glob("nc*_w*_s*")):
        parsed = parse_tag(d.name)
        if not parsed:
            continue
        nc, w, seed = parsed
        rows = load_config(d)
        if rows:
            runs[(nc, w)][seed] = rows

    # Only the (metric, statistic) chosen by the SELECTION stage is evaluated per
    # window. Scanning all nine here and reporting the best would re-introduce exactly
    # the selection-on-evaluation that STUDY.md 3 exists to prevent.
    selected = {}
    sel_path = args.out / "selection_results.json"
    if sel_path.exists():
        sel = json.loads(sel_path.read_text()).get("chosen_per_window") or {}
        for w, r in sel.items():
            selected[int(w)] = (r["metric"], r["statistic"])
    print(f"selected per window (from selection stage): {selected}\n")

    results = []
    for (nc, w), by_seed in sorted(runs.items()):
        seeds = sorted(by_seed)
        if len(seeds) < 2:
            continue
        combos = ([selected[w]] if w in selected
                  else [(m, st) for m in METRICS for st in STATISTICS])
        for metric, stat in combos:
            if True:
                a = collapse(by_seed[seeds[0]], metric, stat)
                b = collapse(by_seed[seeds[1]], metric, stat)
                shared = sorted(set(a) & set(b))
                if len(shared) < 10:
                    continue

                # ---- S1 reproducibility
                rho = spearman([a[k] for k in shared], [b[k] for k in shared])
                deltas = sorted(abs(a[k] - b[k]) for k in shared)
                vals = sorted(a[k] for k in shared)
                spread = vals[int(0.9 * len(vals))] - vals[int(0.1 * len(vals))]

                def sc(ids, table=a):
                    return [table[i] for i in ids if i in table]

                # ---- S2 classification
                s2 = auroc(sc(sets["active"]), sc(sets["scramble"]))

                # ---- E1 potency pair
                hi_s, lo_s = sc(fam.get("YGGFL", [])), sc(fam.get("YGPFL", []))
                e1 = auroc(hi_s, lo_s)
                e1lo, e1hi = bootstrap_auroc(hi_s, lo_s)
                # A windowed screen collapses each family to its distinct windows, so
                # the number of genuinely independent comparisons is the product of the
                # UNIQUE scores, not of the peptide counts. At window 5 both families
                # collapse to one point each: AUROC is then 0 or 1 with no power, and a
                # bootstrap of 1-vs-1 returns a degenerate [1,1] interval. Recorded so
                # that figure can never be read as a result.
                e1_eff = len(set(hi_s)) * len(set(lo_s))
                # ---- E2 potency ranking
                pairs = [(a[i], float(bench[i]["pKi"]))
                         for i in sets["active"] if i in a and bench[i]["pKi"]]
                e2 = spearman([p[0] for p in pairs], [p[1] for p in pairs])
                # ---- E3 added value
                e3 = auroc(sc(sets["active"]), sc(sets["tyr_hard"]))
                e3lo, e3hi = bootstrap_auroc(sc(sets["active"]), sc(sets["tyr_hard"]))
                # ---- E4 held-out ligand retrieval (endomorphin-1 = CHEMBL316446)
                endo = next((i for i in sets["active"] if i.endswith("CHEMBL316446")), None)
                e4 = (percentile_of(a[endo], sc(sets["tyr_hard"]))
                      if endo and endo in a else float("nan"))
                # ---- E5 positional rule
                e5 = auroc(sc(sets["positional_start"]), sc(sets["positional_span"]))

                results.append({
                    "n_conformers": nc, "window": w, "metric": metric,
                    "statistic": stat,
                    "S1_spearman": rho,
                    "S1_median_abs_delta": deltas[len(deltas) // 2],
                    "S1_p10_p90_spread": spread,
                    "S1_pass": bool(rho == rho and rho >= S1_GATE),
                    "S2_auroc_vs_scramble": s2,
                    "E1_potency_pair_auroc": e1, "E1_ci95": [e1lo, e1hi],
                    "E1_unique_high": len(set(hi_s)), "E1_unique_low": len(set(lo_s)),
                    "E1_effective_comparisons": e1_eff,
                    "E1_powered": bool(e1_eff >= 10),
                    "E2_spearman_vs_pKi": e2,
                    "E3_auroc_vs_tyr_hard": e3, "E3_ci95": [e3lo, e3hi],
                    "E4_endomorphin_percentile": e4,
                    "E5_positional_auroc": e5,
                })

    # ---- pre-registered selection: S1 gate, then best S2, ties by cost -----------
    def choose_configuration(rows):
        eligible = [r for r in rows if r["S1_pass"]]
        if not eligible:
            return None
        return max(eligible, key=lambda r: (r["S2_auroc_vs_scramble"],
                                            -r["n_conformers"]))

    chosen = choose_configuration(
        [{k: v for k, v in r.items()
          if k.startswith(("S1", "S2", "n_conf", "window", "metric", "statistic"))}
         for r in results])
    if chosen:
        chosen = next(r for r in results
                      if all(r[k] == chosen[k] for k in
                             ("n_conformers", "window", "metric", "statistic")))

    payload = {"generated": date.today().isoformat(), "s1_gate": S1_GATE,
               "n_configurations": len(results), "chosen": chosen,
               "results": results}
    (args.out / "sweep_results.json").write_text(json.dumps(payload, indent=2) + "\n")

    print(f"{len(results)} configurations evaluated; S1 gate rho >= {S1_GATE}")
    passing = [r for r in results if r["S1_pass"]]
    print(f"passing S1: {len(passing)}")
    hdr = ("  nconf win metric      stat    S1rho   S2auc   E1pair    E2rho   E3auc  E5pos   (! = E1 underpowered)")
    print(hdr)
    for r in sorted(results, key=lambda r: -r["S2_auroc_vs_scramble"])[:12]:
        print("  %5d %3d %-11s %-6s %6.3f %7.3f %7.3f%-3s %6.3f %7.3f %6.3f %s" % (
            r["n_conformers"], r["window"], r["metric"], r["statistic"],
            r["S1_spearman"], r["S2_auroc_vs_scramble"], r["E1_potency_pair_auroc"],
            "" if r["E1_powered"] else "!", r["E2_spearman_vs_pKi"],
            r["E3_auroc_vs_tyr_hard"], r["E5_positional_auroc"],
            "S1ok" if r["S1_pass"] else ""))
    if chosen:
        print("\nCHOSEN (on S1 gate + S2 only):")
        print(json.dumps(chosen, indent=2))
    else:
        print("\nNO configuration passed the S1 reproducibility gate.")


if __name__ == "__main__":
    main()
