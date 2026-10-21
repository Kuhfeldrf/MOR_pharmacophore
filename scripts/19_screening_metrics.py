#!/usr/bin/env python3
"""Early-recognition metrics, and E3 with the confounding property matched.

Both changes come from the methodology literature rather than from this project's own
measurements, and both say the study's existing numbers answer the wrong question.

1. AUROC is the wrong metric for a screen
-----------------------------------------
Truchon & Bayly (J Chem Inf Model 2007, 47:488-508) show that ROC AUC is equivalent to
an average of the ranks of the actives, so good early recognition is cancelled out by
late recognition. For a screen whose entire purpose is the top of the list, that is the
wrong summary. They propose BEDROC, an exponentially weighted generalisation of ROC,
with alpha = 20 so that 80% of the weight falls in the top 8% of the list.

AUROC is NOT wrong everywhere: E1 (YGGFL vs YGPFL) is a two-sample comparison, where
AUROC is the Mann-Whitney U statistic and is exactly the right choice. The metric error
is confined to the screening-shaped contrasts - E3 and the library ranking - where early
recognition is the whole question. Both are computed here.

A caveat from the same paper: BEDROC's standard deviation is only well controlled with
around 50 actives. This study has 36, and 24 reference ligands in the library ranking,
so a permutation null is computed rather than relying on an asymptotic expectation.

2. The decoys are not matched on the property that drives the score
-------------------------------------------------------------------
The `tyr_hard` decoys are matched to the actives on Tyr-start and length, but not on
aromatic count - the property that correlates with the score at rho +0.909. In the
benchmark-design literature that is `artificial enrichment`: performance attributable to
the difference in chemical space between actives and decoys rather than to binding
(Bioinformatics 2021, 37:2134; PLOS ONE 2019, e0220113 on DUD-E's hidden bias). DUD-E's
own principle is property similarity with topological dissimilarity, and this benchmark
satisfies it for length and the N-terminal residue while violating it for exactly the
property that matters here.

So E3 is recomputed with aromatic count matched: within each aromatic-count stratum,
actives against Tyr-starting library decoys of the same stratum. This is the test the
post-hoc stratification of script 14 was reaching for and got wrong, done as a matched
benchmark instead of a score transform. If E3 survives it, the feature-composition
filter of docs/optimisation_study.md 8.4 is real. If it collapses, the screen reduces to
counting aromatic groups and E3's 9x margin was artificial enrichment.

Everything here re-reads scores already on disk: the library screen covered all 14,571
library peptides and the reference ligands at window 8 / 800 conformers / seed 61453,
and the benchmark run used the same configuration and seed, so the two are directly
comparable and nothing needs rescreening.

Usage:
    python scripts/19_screening_metrics.py --library results/library_w8 \\
        --bench-run results/stageB/nc800_w8_s61453 --out docs
"""
import argparse
import csv
import glob
import json
import random
from collections import defaultdict
from datetime import date
from math import exp, sinh, cosh
from pathlib import Path

METRIC = "tversky_ref_mean"
ALPHA = 20.0
AROMATIC = "FWYH"
EF_FRACTIONS = (0.01, 0.05, 0.10)


def n_aromatic(w):
    return sum(w.count(c) for c in AROMATIC)


# ---------------------------------------------------------------- metrics
def auroc(pos, neg):
    if not pos or not neg:
        return float("nan")
    wins = sum(1.0 if p > n else (0.5 if p == n else 0.0) for p in pos for n in neg)
    return wins / (len(pos) * len(neg))


def rie_bedroc(ranks, n_total, n_actives, alpha=ALPHA):
    """RIE and BEDROC from 1-based ranks of the actives (Truchon & Bayly 2007)."""
    N, n = float(n_total), float(n_actives)
    if not n or N <= 1:
        return float("nan"), float("nan")
    s = sum(exp(-alpha * r / N) for r in ranks)
    denom = (1.0 / N) * (1.0 - exp(-alpha)) / (exp(alpha / N) - 1.0)
    rie = (s / n) / denom
    ra = n / N
    scale = ra * sinh(alpha / 2.0) / (cosh(alpha / 2.0) - cosh(alpha / 2.0 - alpha * ra))
    bedroc = rie * scale + 1.0 / (1.0 - exp(alpha * (1.0 - ra)))
    return rie, bedroc


def enrichment(ranks, n_total, n_actives, frac):
    cut = max(1, int(round(frac * n_total)))
    hits = sum(1 for r in ranks if r <= cut)
    return (hits / cut) / (n_actives / n_total) if n_actives else float("nan")


def ranks_of(active_ids, scores):
    """1-based ranks of the actives in a descending score ordering, ties averaged."""
    order = sorted(scores, key=lambda k: -scores[k])
    rank = {}
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and scores[order[j + 1]] == scores[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            rank[order[k]] = avg
        i = j + 1
    return [rank[a] for a in active_ids if a in rank]


def permutation_null(n_total, n_actives, n_iter=2000, seed=11, alpha=ALPHA):
    """Empirical null for BEDROC, since 36 actives is below the ~50 recommended."""
    rng = random.Random(seed)
    vals = []
    for _ in range(n_iter):
        rs = rng.sample(range(1, n_total + 1), n_actives)
        vals.append(rie_bedroc(rs, n_total, n_actives, alpha)[1])
    vals.sort()
    return {"median": vals[len(vals) // 2],
            "p95": vals[int(0.95 * len(vals))],
            "p99": vals[int(0.99 * len(vals))]}


# ---------------------------------------------------------------- loading
def load_blocks(d):
    rows = {}
    for f in sorted(glob.glob(str(Path(d) / "block*.tsv"))):
        with open(f) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                rows[r["id"]] = r
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--library", default=Path("results/library_w8"), type=Path)
    ap.add_argument("--bench-run", default=Path("results/stageB/nc800_w8_s61453"),
                    type=Path)
    ap.add_argument("--benchmark", default=Path("data/benchmark/benchmark.tsv"),
                    type=Path)
    ap.add_argument("--out", default=Path("docs"), type=Path)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    lib_rows = load_blocks(args.library)
    bench_rows = load_blocks(args.bench_run)
    with args.benchmark.open() as fh:
        bench = {r["id"]: r for r in csv.DictReader(fh, delimiter="\t")}

    out = {"generated": date.today().isoformat(), "metric": METRIC, "alpha": ALPHA,
           "configuration": {"window": 8, "n_conformers": 800, "seed": 61453}}

    # ================= A. retrieval of the known ligands from the library ===========
    # The real screening question: ranking all 14,571 library peptides plus the known
    # reference ligands by score, how early do the known ligands appear?
    scores, refs = {}, []
    for pid, r in lib_rows.items():
        scores[pid] = float(r[METRIC])
        if r["set"] != "library":
            refs.append(pid)
    endogenous = [p for p in refs
                  if lib_rows[p]["set"] in ("positive_endogenous", "template_ligand")]

    sect = {}
    for label, actives in (("all_reference_ligands", refs),
                           ("endogenous_and_template_only", endogenous)):
        rk = ranks_of(actives, scores)
        N, n = len(scores), len(rk)
        rie, bed = rie_bedroc(rk, N, n)
        null = permutation_null(N, n)
        sect[label] = {
            "n_actives": n, "n_total": N,
            "auroc": auroc([scores[a] for a in actives if a in scores],
                           [v for k, v in scores.items() if k not in set(actives)]),
            "RIE": rie, "BEDROC_alpha20": bed,
            "BEDROC_permutation_null": null,
            "BEDROC_exceeds_null_p99": bool(bed > null["p99"]),
            "enrichment_factor": {str(f): enrichment(rk, N, n, f)
                                  for f in EF_FRACTIONS},
            "median_rank": sorted(rk)[len(rk) // 2],
            "best_rank": min(rk), "worst_rank": max(rk),
        }
    out["A_library_retrieval"] = sect

    # ================= B. E3 recomputed with aromatic count MATCHED ================
    # Actives from the benchmark run; decoys are Tyr-starting library peptides of the
    # same window aromatic count. Distinct windows only, so one repeat family cannot
    # dominate a stratum.
    act = {}
    for bid, m in bench.items():
        if m["set"] == "active" and bid in bench_rows:
            r = bench_rows[bid]
            act[bid] = (float(r[METRIC]), n_aromatic(r["window"]))
    decoy_by_arom = defaultdict(dict)
    for pid, r in lib_rows.items():
        if r["set"] != "library" or not r["sequence"].startswith("Y"):
            continue
        decoy_by_arom[n_aromatic(r["window"])].setdefault(
            r["window"], float(r[METRIC]))

    strata, pooled_pos, pooled_neg = [], [], []
    num = den = 0.0
    for a in sorted({v[1] for v in act.values()}):
        pos = [s for s, ar in act.values() if ar == a]
        neg = list(decoy_by_arom.get(a, {}).values())
        if not pos or len(neg) < 5:
            strata.append({"aromatic_count": a, "n_actives": len(pos),
                           "n_decoys": len(neg), "auroc": None,
                           "skipped": "fewer than 5 matched decoy windows"})
            continue
        u = auroc(pos, neg)
        w = len(pos) * len(neg)
        num += u * w
        den += w
        pooled_pos += pos
        pooled_neg += neg
        strata.append({"aromatic_count": a, "n_actives": len(pos),
                       "n_decoys": len(neg), "auroc": u,
                       "median_active": sorted(pos)[len(pos) // 2],
                       "median_decoy": sorted(neg)[len(neg) // 2]})
    out["B_E3_aromatic_matched"] = {
        "per_stratum": strata,
        "pooled_auroc_weighted": num / den if den else float("nan"),
        "unmatched_E3_for_comparison": 0.733,
        "note": ("Pooling is weighted by stratum comparison count. The unmatched E3 of "
                 "0.733 is the figure this replaces; the difference between them is "
                 "the artificial enrichment attributable to aromatic count."),
    }

    # ================= C. BEDROC for E3, matched and unmatched =====================
    def bedroc_for(pos_vals, neg_vals):
        tbl = {}
        for i, v in enumerate(pos_vals):
            tbl[f"a{i}"] = v
        for i, v in enumerate(neg_vals):
            tbl[f"d{i}"] = v
        rk = ranks_of([f"a{i}" for i in range(len(pos_vals))], tbl)
        N, n = len(tbl), len(rk)
        rie, bed = rie_bedroc(rk, N, n)
        return {"RIE": rie, "BEDROC_alpha20": bed, "n_total": N, "n_actives": n,
                "BEDROC_permutation_null": permutation_null(N, n),
                "enrichment_factor": {str(f): enrichment(rk, N, n, f)
                                      for f in EF_FRACTIONS}}

    unmatched_neg = [float(bench_rows[b][METRIC]) for b, m in bench.items()
                     if m["set"] == "tyr_hard" and b in bench_rows]
    out["C_E3_bedroc"] = {
        "unmatched_decoys": bedroc_for([v[0] for v in act.values()], unmatched_neg),
        "aromatic_matched_decoys": bedroc_for(pooled_pos, pooled_neg),
    }

    (args.out / "screening_metrics.json").write_text(json.dumps(out, indent=2) + "\n")

    # ---------------------------------------------------------------- report
    print("A. Retrieval of known ligands from the full library ranking")
    print("   (the actual screening question; BEDROC alpha=20)\n")
    for label, d in sect.items():
        print("  %s: %d actives in %d ranked" % (label, d["n_actives"], d["n_total"]))
        print("     AUROC  %.3f" % d["auroc"])
        print("     BEDROC %.3f   permutation null median %.3f, p99 %.3f  -> %s"
              % (d["BEDROC_alpha20"], d["BEDROC_permutation_null"]["median"],
                 d["BEDROC_permutation_null"]["p99"],
                 "above null" if d["BEDROC_exceeds_null_p99"] else "NOT above null"))
        print("     RIE    %.2f" % d["RIE"])
        print("     EF     " + "  ".join(
            "%s%%: %.1f" % (int(float(f) * 100), v)
            for f, v in d["enrichment_factor"].items()))
        print("     ranks  best %d, median %d, worst %d\n"
              % (d["best_rank"], d["median_rank"], d["worst_rank"]))

    print("B. E3 with aromatic count MATCHED between actives and decoys")
    print("   arom  actives  decoys   AUROC   median act / decoy")
    for s in strata:
        if s.get("auroc") is None:
            print("   %4d  %7d  %6d   %s" % (s["aromatic_count"], s["n_actives"],
                                             s["n_decoys"], s["skipped"]))
        else:
            print("   %4d  %7d  %6d   %.3f   %.4f / %.4f"
                  % (s["aromatic_count"], s["n_actives"], s["n_decoys"], s["auroc"],
                     s["median_active"], s["median_decoy"]))
    print("\n   pooled (comparison-weighted) : %.3f"
          % out["B_E3_aromatic_matched"]["pooled_auroc_weighted"])
    print("   unmatched E3 it replaces     : 0.733")

    print("\nC. E3 as BEDROC, both decoy sets")
    for k, d in out["C_E3_bedroc"].items():
        print("   %-26s BEDROC %.3f (null median %.3f, p99 %.3f)  EF1%% %.1f"
              % (k, d["BEDROC_alpha20"], d["BEDROC_permutation_null"]["median"],
                 d["BEDROC_permutation_null"]["p99"],
                 d["enrichment_factor"]["0.01"]))
    print("\nwrote %s" % (args.out / "screening_metrics.json"))


if __name__ == "__main__":
    main()
