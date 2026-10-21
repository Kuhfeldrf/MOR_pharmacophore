#!/usr/bin/env python3
"""Rank and tier the full library screen at the validated configuration.

Tiers are defined from the NEGATIVE distribution (the Tyr-matched decoys of the
benchmark), not from the library's own deciles. The first study showed why: library
deciles describe the conformer seed as much as the peptide, because the library score
distribution is compressed relative to the seed jitter. Anchoring the boundaries to a
background whose members are presumed inactive keeps a tier meaning the same thing when
the method is pointed at another receptor.

Every boundary is checked against the measured seed jitter. A tier edge closer to its
neighbour than the jitter is reported as unsafe rather than used.

The free-N-terminus rule is reported as a separate flag and is never folded into the
score, so its contribution stays separable from the 3D matching's.

Usage:
    python scripts/13_library_report.py --library results/library_w8 \\
        --benchmark-run results/stageB/nc800_w8_s61453 --out docs
"""
import argparse
import csv
import json
from collections import Counter
from datetime import date
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from mor_pharmacophore.alignit import read_scores  # noqa: E402

METRIC = "tversky_ref"
STATISTIC = "mean"
# Tier edges as percentiles of the decoy background.
TIER_EDGES = [("tier1", 0.99), ("tier2", 0.95), ("tier3", 0.90)]


def collapse(rows, metric=METRIC, statistic=STATISTIC):
    grouped = {}
    for r in rows:
        grouped.setdefault(r["candidate"], []).append(r[metric])
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


def load_run(d):
    rows = []
    for tab in sorted(Path(d).glob("block*/scores.tab")):
        rows.extend(read_scores(tab))
    return rows


def percentile(sorted_vals, q):
    if not sorted_vals:
        return float("nan")
    i = min(len(sorted_vals) - 1, max(0, int(q * len(sorted_vals))))
    return sorted_vals[i]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--library", default=Path("results/library_w8"), type=Path)
    ap.add_argument("--benchmark-run",
                    default=Path("results/stageB/nc800_w8_s61453"), type=Path)
    ap.add_argument("--benchmark", default=Path("data/benchmark/benchmark.tsv"),
                    type=Path)
    ap.add_argument("--targets", default=Path("data/targets.tsv"), type=Path)
    ap.add_argument("--library-tsv", default=Path("data/library/peptides.tsv"),
                    type=Path)
    ap.add_argument("--jitter", default=None, type=float,
                    help="measured median |delta| between seeds; read from "
                         "sweep_results.json when omitted")
    ap.add_argument("--out", default=Path("docs"), type=Path)
    ap.add_argument("--top", default=50, type=int)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    # ---- the decoy background, scored under the same configuration ---------------
    with args.benchmark.open() as fh:
        bench = {r["id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    bench_scores = collapse(load_run(args.benchmark_run))
    decoys = sorted(v for k, v in bench_scores.items()
                    if bench.get(k, {}).get("set") == "tyr_hard")
    actives = sorted(v for k, v in bench_scores.items()
                     if bench.get(k, {}).get("set") == "active")
    if not decoys:
        raise SystemExit("no decoy scores found; cannot define tiers")

    jitter = args.jitter
    if jitter is None:
        sweep = args.out / "sweep_results.json"
        if sweep.exists():
            for r in json.loads(sweep.read_text())["results"]:
                if r["window"] == 8 and r["n_conformers"] == 800:
                    jitter = r["S1_median_abs_delta"]
        if jitter is None:
            raise SystemExit("--jitter not given and not found in sweep_results.json")

    edges = [(name, percentile(decoys, q), q) for name, q in TIER_EDGES]
    unsafe = [edges[i][0] for i in range(len(edges) - 1)
              if abs(edges[i][1] - edges[i + 1][1]) < jitter]

    # ---- library scores -----------------------------------------------------------
    lib_scores = collapse(load_run(args.library))
    with args.targets.open() as fh:
        targets = {r["id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    with args.library_tsv.open() as fh:
        meta = {r["peptide_id"]: r for r in csv.DictReader(fh, delimiter="\t")}

    def tier_of(score):
        for name, cut, _ in edges:
            if score >= cut:
                return name
        return "untiered"

    rows = []
    for pid, score in lib_scores.items():
        t = targets.get(pid, {})
        if t.get("set") != "library":
            continue
        m = meta.get(pid, {})
        seq = t.get("sequence", "")
        rows.append({
            "peptide_id": pid, "sequence": seq, "length": len(seq),
            "score": score, "tier": tier_of(score),
            "free_nterm_tyr": "yes" if seq.startswith("Y") else "no",
            "phases": m.get("phases", ""), "matrices": m.get("matrices", ""),
        })
    rows.sort(key=lambda r: -r["score"])

    out_tsv = args.out / "library_ranked.tsv"
    cols = ["rank", "peptide_id", "sequence", "length", "score", "tier",
            "free_nterm_tyr", "phases", "matrices"]
    with out_tsv.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for i, r in enumerate(rows, 1):
            fh.write("\t".join(str(x) for x in [
                i, r["peptide_id"], r["sequence"], r["length"],
                f"{r['score']:.4f}", r["tier"], r["free_nterm_tyr"],
                r["phases"], r["matrices"]]) + "\n")

    tiers = Counter(r["tier"] for r in rows)
    tyr_by_tier = {t: sum(1 for r in rows
                          if r["tier"] == t and r["free_nterm_tyr"] == "yes")
                   for t in tiers}
    summary = {
        "generated": date.today().isoformat(),
        "configuration": {"window": 8, "n_conformers": 800,
                          "metric": METRIC, "statistic": STATISTIC},
        "n_library_scored": len(rows),
        "seed_jitter_median_abs_delta": jitter,
        "tier_edges": [{"tier": n, "decoy_percentile": q, "score_cut": c}
                       for n, c, q in edges],
        "unsafe_tier_edges": unsafe,
        "tier_counts": dict(tiers),
        "tyr_starting_by_tier": tyr_by_tier,
        "decoy_background_n": len(decoys),
        "active_median": actives[len(actives) // 2] if actives else None,
        "decoy_median": decoys[len(decoys) // 2],
        "top": [{k: r[k] for k in ("peptide_id", "sequence", "score", "tier",
                                   "free_nterm_tyr", "phases")}
                for r in rows[:args.top]],
    }
    (args.out / "library_report.json").write_text(json.dumps(summary, indent=2) + "\n")

    print(f"scored {len(rows)} library peptides at window 8 / 800 conformers")
    print(f"seed jitter (median |delta|): {jitter:.4f}")
    for n, c, q in edges:
        print(f"  {n}: score >= {c:.4f}  (decoy p{int(q*100)})")
    if unsafe:
        print(f"  WARNING unsafe tier edges (closer than the jitter): {unsafe}")
    for t in ("tier1", "tier2", "tier3", "untiered"):
        if t in tiers:
            print(f"  {t:9s} {tiers[t]:6d}   Tyr-starting {tyr_by_tier.get(t,0)}")
    print(f"\nwrote {out_tsv}")
    print("\nTop 15:")
    for i, r in enumerate(rows[:15], 1):
        print("  %2d %-14s %-26s %.4f %-8s Y=%s %s" % (
            i, r["peptide_id"], r["sequence"][:26], r["score"], r["tier"],
            r["free_nterm_tyr"], r["phases"]))


if __name__ == "__main__":
    main()
