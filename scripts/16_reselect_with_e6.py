#!/usr/bin/env python3
"""Re-select the scoring metric with E6 (composition confound) as a validity gate.

Why this exists
---------------
The original selection stage (`11_select.py`) chose TVERSKY_REF on S1 and S2 alone.
Both are blind to the confound found at library scale: S2's negatives are
composition-matched scrambles and S1 is a seed-agreement check, so neither can see that
the score tracks aromatic-residue count (rho +0.909 on the library under TVERSKY_REF).

TVERSKY_REF is asymmetric - it measures how much of the REFERENCE pharmacophore the
candidate covers and never penalises features the candidate has in excess. A His-rich
window therefore accumulates coverage (His contributes ring + donor + acceptor at one
residue) without presenting the message motif. TANIMOTO is symmetric and TVERSKY_DB is
asymmetric the other way, so both should penalise excess features. That is a prediction
about the mechanism, and it is testable on score tables already on disk.

Status of each number here
--------------------------
  S1  selection metric, unchanged, gate rho >= 0.95
  S2  selection metric, unchanged
  E6  VALIDITY GATE, newly pre-registered in STUDY.md 7: Spearman rho between score and
      aromatic-residue count over the screened library, ceiling 0.6. A gate is a
      statement about whether a measurement means anything, not about how well the
      method performs, which is why using it to exclude a metric is not
      selection-on-evaluation. It is computed on the LIBRARY, which contains no
      benchmark peptide and no affinity data, so it cannot encode E1-E5.
  E1-E5  evaluation only, reported for every surviving combination, never used to choose.

The honest risk, stated up front: E6 was defined after seeing that TVERSKY_REF failed
it. Choosing a different metric on a gate defined after the fact is weaker evidence than
pre-registration, so E1-E5 are reported for ALL nine combinations rather than only the
winner, and the full table is written out. If a metric that passes E6 also happens to
hold E1 and E3, that is a real improvement; if it buys E6 by destroying E1, that is the
finding and TVERSKY_REF stays with the post-hoc shortlist of script 15.

Usage:
    python scripts/16_reselect_with_e6.py --library results/library_w8 \\
        --bench-a results/stageB/nc800_w8_s20261020 \\
        --bench-b results/stageB/nc800_w8_s61453 --out docs
"""
import argparse
import csv
import glob
import json
import random
from collections import defaultdict
from datetime import date
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from mor_pharmacophore.alignit import read_scores  # noqa: E402

METRICS = ("tanimoto", "tversky_ref", "tversky_db")
STATISTICS = ("max", "mean", "top10")
S1_GATE = 0.95
E6_CEILING = 0.6
AROMATIC = "FWYH"


# ---------------------------------------------------------------- statistics
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


def auroc(pos, neg):
    if not pos or not neg:
        return float("nan")
    wins = sum(1.0 if p > n else (0.5 if p == n else 0.0) for p in pos for n in neg)
    return wins / (len(pos) * len(neg))


def bootstrap_auroc(pos, neg, n=2000, seed=7):
    if not pos or not neg:
        return float("nan"), float("nan")
    rng = random.Random(seed)
    vals = sorted(auroc([rng.choice(pos) for _ in pos],
                        [rng.choice(neg) for _ in neg]) for _ in range(n))
    return vals[int(0.025 * n)], vals[int(0.975 * n) - 1]


def n_aromatic(window):
    return sum(window.count(c) for c in AROMATIC)


# ---------------------------------------------------------------- loading
def collapse_all(run_dir):
    """candidate -> metric -> statistic, in one pass over the per-conformer tables.

    The library tables total ~1.2 GB, so the per-conformer values are accumulated
    rather than held: sum and count for the mean, a bounded top-k heap-free list for
    max and the top decile. Keeping the best 10% of a 800-conformer ensemble is 80
    floats per candidate per metric, which is affordable; keeping all of them is not.
    """
    acc = defaultdict(lambda: {m: {"sum": 0.0, "n": 0, "top": []} for m in METRICS})
    for tab in sorted(glob.glob(str(Path(run_dir) / "block*" / "scores.tab"))):
        for r in read_scores(tab):
            a = acc[r["candidate"]]
            for m in METRICS:
                v = r[m]
                d = a[m]
                d["sum"] += v
                d["n"] += 1
                t = d["top"]
                t.append(v)
                if len(t) > 256:                 # keep only the running best 128
                    t.sort(reverse=True)
                    del t[128:]
    out = {}
    for cand, a in acc.items():
        per = {}
        for m in METRICS:
            d = a[m]
            if not d["n"]:
                continue
            t = sorted(d["top"], reverse=True)
            k = max(1, d["n"] // 10)
            per[m] = {
                "max": t[0],
                "mean": d["sum"] / d["n"],
                "top10": sum(t[:k]) / min(k, len(t)),
            }
        out[cand] = per
    return out


def windows_of(run_dir):
    """candidate id -> scored window, from the per-peptide block tables."""
    w = {}
    for f in sorted(glob.glob(str(Path(run_dir) / "block*.tsv"))):
        with open(f) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                w[r["id"]] = (r["window"], r["set"], r["sequence"])
    return w


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--library", default=Path("results/library_w8"), type=Path)
    ap.add_argument("--bench-a", default=Path("results/stageB/nc800_w8_s20261020"),
                    type=Path)
    ap.add_argument("--bench-b", default=Path("results/stageB/nc800_w8_s61453"),
                    type=Path)
    ap.add_argument("--benchmark", default=Path("data/benchmark/benchmark.tsv"),
                    type=Path)
    ap.add_argument("--out", default=Path("docs"), type=Path)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    print("loading library score tables ...", flush=True)
    lib = collapse_all(args.library)
    lib_meta = windows_of(args.library)
    print("loading benchmark score tables (two seeds) ...", flush=True)
    ba, bb = collapse_all(args.bench_a), collapse_all(args.bench_b)

    with args.benchmark.open() as fh:
        bench = {r["id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    sets = defaultdict(list)
    fam = defaultdict(list)
    for bid, r in bench.items():
        sets[r["set"]].append(bid)
        if r["set"] == "active" and r["family"]:
            fam[r["family"]].append(bid)

    # distinct library windows, for E6 (one point per window, so a repeat family that
    # contributes 30 peptides with one window is not counted 30 times)
    lib_windows = {}
    for cand, per in lib.items():
        meta = lib_meta.get(cand)
        if not meta or meta[1] != "library":
            continue
        lib_windows.setdefault(meta[0], per)

    results = []
    for metric in METRICS:
        for stat in STATISTICS:
            # ---- E6 validity gate, on the library
            ws = sorted(lib_windows)
            sc = [lib_windows[w][metric][stat] for w in ws]
            e6 = spearman(sc, [n_aromatic(w) for w in ws])
            e6_his = spearman(sc, [w.count("H") for w in ws])

            # ---- S1 / S2 on the benchmark
            a = {c: p[metric][stat] for c, p in ba.items() if metric in p}
            b = {c: p[metric][stat] for c, p in bb.items() if metric in p}
            shared = sorted(set(a) & set(b))
            rho = spearman([a[k] for k in shared], [b[k] for k in shared])
            deltas = sorted(abs(a[k] - b[k]) for k in shared)
            vals = sorted(a[k] for k in shared)
            spread = vals[int(0.9 * len(vals))] - vals[int(0.1 * len(vals))]

            def g(ids, table=b):
                return [table[i] for i in ids if i in table]

            s2 = auroc(g(sets["active"]), g(sets["scramble"]))
            hi, lo = g(fam.get("YGGFL", [])), g(fam.get("YGPFL", []))
            act, hard = g(sets["active"]), g(sets["tyr_hard"])
            e1lo, e1hi = bootstrap_auroc(hi, lo)
            e3lo, e3hi = bootstrap_auroc(act, hard)
            endo = next((i for i in sets["active"] if i.endswith("CHEMBL316446")), None)
            pairs = [(b[i], float(bench[i]["pKi"]))
                     for i in sets["active"] if i in b and bench[i]["pKi"]]

            results.append({
                "metric": metric, "statistic": stat,
                "S1_spearman": rho,
                "S1_median_abs_delta": deltas[len(deltas) // 2],
                "S1_pass": bool(rho == rho and rho >= S1_GATE),
                "S2_auroc_vs_scramble": s2,
                "E6_rho_aromatic": e6, "E6_rho_his": e6_his,
                "E6_pass": bool(e6 == e6 and abs(e6) <= E6_CEILING),
                "E1_potency_pair_auroc": auroc(hi, lo), "E1_ci95": [e1lo, e1hi],
                "E1_effective_comparisons": len(set(hi)) * len(set(lo)),
                "E2_spearman_vs_pKi": spearman([p[0] for p in pairs],
                                               [p[1] for p in pairs]),
                "E3_auroc_vs_tyr_hard": auroc(act, hard), "E3_ci95": [e3lo, e3hi],
                "E4_endomorphin_percentile": (
                    sum(1 for v in hard if v <= b[endo]) / len(hard)
                    if endo in b and hard else float("nan")),
                "E5_positional_auroc": auroc(g(sets["positional_start"]),
                                             g(sets["positional_span"])),
                "n_library_windows": len(ws), "S1_p10_p90_spread": spread,
            })

    eligible = [r for r in results if r["S1_pass"] and r["E6_pass"]]
    chosen = max(eligible, key=lambda r: r["S2_auroc_vs_scramble"]) if eligible else None

    payload = {
        "generated": date.today().isoformat(),
        "configuration": {"window": 8, "n_conformers": 800},
        "gates": {"S1_spearman_min": S1_GATE, "E6_abs_rho_max": E6_CEILING},
        "note": ("E6 is a validity gate computed on the library, which holds no "
                 "benchmark peptide and no affinity data. E1-E5 are reported for all "
                 "nine combinations so the choice can be audited, and were not used "
                 "to choose."),
        "chosen": chosen, "results": results,
    }
    (args.out / "reselection_e6.json").write_text(json.dumps(payload, indent=2) + "\n")

    print("\nAll nine metric x statistic combinations at window 8 / 800 conformers")
    print("gates: S1 rho >= %.2f, |E6 rho| <= %.2f\n" % (S1_GATE, E6_CEILING))
    print("  metric      stat    S1rho  S2auc |  E6arom E6his | E1pair         E2rho "
          " E3auc          E4    E5     gates")
    for r in sorted(results, key=lambda r: (r["metric"], r["statistic"])):
        print("  %-11s %-6s %6.3f %6.3f | %+6.3f %+6.3f | %5.3f [%.2f,%.2f] %6.3f "
              "%5.3f [%.2f,%.2f] %5.3f %5.3f  %s%s" % (
                  r["metric"], r["statistic"], r["S1_spearman"],
                  r["S2_auroc_vs_scramble"], r["E6_rho_aromatic"], r["E6_rho_his"],
                  r["E1_potency_pair_auroc"], r["E1_ci95"][0], r["E1_ci95"][1],
                  r["E2_spearman_vs_pKi"], r["E3_auroc_vs_tyr_hard"],
                  r["E3_ci95"][0], r["E3_ci95"][1],
                  r["E4_endomorphin_percentile"], r["E5_positional_auroc"],
                  "S1 " if r["S1_pass"] else "-- ", "E6" if r["E6_pass"] else "--"))
    if chosen:
        print("\nPASSES BOTH GATES, highest S2: %s / %s" % (chosen["metric"],
                                                            chosen["statistic"]))
        print("  E6 rho vs aromatic count : %+.3f" % chosen["E6_rho_aromatic"])
        print("  E1 potency pair          : %.3f %s" % (
            chosen["E1_potency_pair_auroc"],
            [round(x, 3) for x in chosen["E1_ci95"]]))
        print("  E3 vs Tyr-matched decoys : %.3f %s" % (
            chosen["E3_auroc_vs_tyr_hard"],
            [round(x, 3) for x in chosen["E3_ci95"]]))
    else:
        print("\nNO combination passes both gates. TVERSKY_REF plus the within-stratum "
              "shortlist (script 15) remains the reported result.")


if __name__ == "__main__":
    main()
