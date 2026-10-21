#!/usr/bin/env python3
"""Post-hoc diagnostic: how much of the library ranking is composition, not shape?

This is NOT one of STUDY.md's pre-registered metrics. It was written after the full
library screen, because the library ranking itself raised the question: the top of the
list is a single His-rich repeat family and no known active reaches tier 1, while the
pre-registered benchmark metrics (E1-E5) all looked reasonable. Labelled post-hoc and
reported as such.

The diagnostic: Spearman rho between the score and simple composition counts of the
scored window, over all 14,571 library windows. If aromatic count alone predicts the
score, the screen is ranking composition and the 3D superposition is decoration.

The remedy: a composition-stratified score. Each peptide's score is converted to its
percentile within the library windows that have the SAME aromatic-residue count, so
"how well does this peptide match the template, for a peptide with this many aromatic
groups to match with" replaces "how many aromatic groups does it have". The library is
the background because it is large and presumed inactive.

The test of the remedy is the pre-registered benchmark, rescored: E1 should be
UNCHANGED (the YGGFL/YGPFL families are exactly aromatic-matched, so stratification
cannot touch that contrast) and E3/E4 should improve if the confound was real and the
stratification removes it. E1 being unchanged is the control on the remedy itself.

Usage:
    python scripts/14_composition_control.py \\
        --library results/library_w8 --benchmark-run results/stageB/nc800_w8_s61453 \\
        --out docs
"""
import argparse
import csv
import glob
import json
import random
from collections import defaultdict
from datetime import date
from pathlib import Path

SCORE = "tversky_ref_mean"
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


# ---------------------------------------------------------------- loading
def load_blocks(d, want=None):
    """Per-peptide collapsed rows (block*.tsv), not the per-conformer tables."""
    out = []
    for f in sorted(glob.glob(str(Path(d) / "block*.tsv"))):
        with open(f) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if want is None or r["set"] in want:
                    out.append(r)
    return out


def n_aromatic(window):
    return sum(window.count(c) for c in AROMATIC)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--library", default=Path("results/library_w8"), type=Path)
    ap.add_argument("--benchmark-run",
                    default=Path("results/stageB/nc800_w8_s61453"), type=Path)
    ap.add_argument("--benchmark", default=Path("data/benchmark/benchmark.tsv"),
                    type=Path)
    ap.add_argument("--out", default=Path("docs"), type=Path)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    lib = [r for r in load_blocks(args.library) if r["set"] == "library"]
    refs = [r for r in load_blocks(args.library) if r["set"] != "library"]

    # ---- 1. the diagnostic ------------------------------------------------------
    s = [float(r[SCORE]) for r in lib]
    w = [r["window"] for r in lib]
    confounds = {
        "aromatic_FWYH": [n_aromatic(x) for x in w],
        "his": [x.count("H") for x in w],
        "matched_features": [float(r["matched_features"]) for r in lib],
        "window_length": [len(x) for x in w],
        "charged_KRDEH": [sum(x.count(c) for c in "KRDEH") for x in w],
        "starts_with_Y": [1 if x.startswith("Y") else 0 for x in w],
    }
    diagnostic = {k: spearman(s, v) for k, v in confounds.items()}

    # ---- 2. the composition-stratified score ------------------------------------
    # Background = distinct library windows, binned by aromatic count. Distinct
    # windows, not peptides, so a repeat family that contributes 30 peptides with one
    # window does not weight the background 30 times.
    bg = defaultdict(list)
    seen = set()
    for r in lib:
        if r["window"] in seen:
            continue
        seen.add(r["window"])
        bg[n_aromatic(r["window"])].append(float(r[SCORE]))
    for k in bg:
        bg[k].sort()

    def stratified(window, score):
        pop = bg.get(n_aromatic(window))
        if not pop or len(pop) < 20:          # too thin a stratum to calibrate against
            return None
        return sum(1 for v in pop if v <= score) / len(pop)

    bin_sizes = {k: len(v) for k, v in sorted(bg.items())}

    # ---- 3. does the remedy hold up on the pre-registered benchmark? -------------
    with args.benchmark.open() as fh:
        bench = {r["id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    brows = load_blocks(args.benchmark_run)
    raw, strat, fam = {}, {}, defaultdict(list)
    for r in brows:
        bid = r["id"]
        if bid not in bench:
            continue
        raw[bid] = float(r[SCORE])
        st = stratified(r["window"], float(r[SCORE]))
        if st is not None:
            strat[bid] = st
    sets = defaultdict(list)
    for bid, r in bench.items():
        sets[r["set"]].append(bid)
        if r["set"] == "active" and r["family"]:
            fam[r["family"]].append(bid)

    def metrics(table):
        def sc(ids):
            return [table[i] for i in ids if i in table]
        hi, lo = sc(fam.get("YGGFL", [])), sc(fam.get("YGPFL", []))
        act, hard = sc(sets["active"]), sc(sets["tyr_hard"])
        e1lo, e1hi = bootstrap_auroc(hi, lo)
        e3lo, e3hi = bootstrap_auroc(act, hard)
        endo = next((i for i in sets["active"] if i.endswith("CHEMBL316446")), None)
        pairs = [(table[i], float(bench[i]["pKi"]))
                 for i in sets["active"] if i in table and bench[i]["pKi"]]
        return {
            "E1_potency_pair_auroc": auroc(hi, lo), "E1_ci95": [e1lo, e1hi],
            "E1_effective_comparisons": len(set(hi)) * len(set(lo)),
            "E2_spearman_vs_pKi": spearman([p[0] for p in pairs],
                                           [p[1] for p in pairs]),
            "E3_auroc_vs_tyr_hard": auroc(act, hard), "E3_ci95": [e3lo, e3hi],
            "E4_endomorphin_percentile": (
                sum(1 for v in hard if v <= table[endo]) / len(hard)
                if endo in table and hard else float("nan")),
            "E5_positional_auroc": auroc(sc(sets["positional_start"]),
                                         sc(sets["positional_span"])),
            "n_scored": len(table),
        }

    m_raw, m_strat = metrics(raw), metrics(strat)

    # ---- 4. re-tier the library on the stratified score --------------------------
    lib_strat = []
    for r in lib:
        st = stratified(r["window"], float(r[SCORE]))
        if st is None:
            continue
        lib_strat.append({"peptide_id": r["id"], "sequence": r["sequence"],
                          "window": r["window"], "raw": float(r[SCORE]),
                          "stratified": st, "n_aromatic": n_aromatic(r["window"])})
    lib_strat.sort(key=lambda r: -r["stratified"])

    ref_strat = sorted(
        ({"id": r["id"].replace("ref__", ""), "set": r["set"],
          "window": r["window"], "raw": float(r[SCORE]),
          "stratified": stratified(r["window"], float(r[SCORE]))}
         for r in refs),
        key=lambda r: -(r["stratified"] if r["stratified"] is not None else -1))

    out_tsv = args.out / "library_ranked_stratified.tsv"
    with out_tsv.open("w") as fh:
        fh.write("rank\tpeptide_id\tsequence\twindow\tn_aromatic\t"
                 "raw_score\tstratified\n")
        for i, r in enumerate(lib_strat, 1):
            fh.write("%d\t%s\t%s\t%s\t%d\t%.4f\t%.4f\n" % (
                i, r["peptide_id"], r["sequence"], r["window"],
                r["n_aromatic"], r["raw"], r["stratified"]))

    payload = {
        "generated": date.today().isoformat(),
        "post_hoc": True,
        "note": "Diagnostic and remedy written after the library screen; not a "
                "pre-registered metric. E1 is the control on the remedy: the "
                "YGGFL/YGPFL families are aromatic-matched, so stratification must "
                "leave E1 unchanged.",
        "score_field": SCORE,
        "confound_spearman": diagnostic,
        "aromatic_bin_sizes": bin_sizes,
        "benchmark_raw": m_raw,
        "benchmark_stratified": m_strat,
        "references_stratified": ref_strat,
        "top_stratified": lib_strat[:50],
    }
    (args.out / "composition_control.json").write_text(
        json.dumps(payload, indent=2) + "\n")

    # ---- report -------------------------------------------------------------------
    print("Diagnostic: Spearman rho of %s against composition (n=%d library windows)"
          % (SCORE, len(lib)))
    for k, v in sorted(diagnostic.items(), key=lambda kv: -abs(kv[1])):
        print("  %-20s %+.3f" % (k, v))
    print("\nAromatic-count strata (distinct windows): %s" % bin_sizes)
    print("\nPre-registered benchmark, raw vs composition-stratified score:")
    print("  %-28s %10s %10s" % ("", "raw", "stratified"))
    for k in ("E1_potency_pair_auroc", "E2_spearman_vs_pKi", "E3_auroc_vs_tyr_hard",
              "E4_endomorphin_percentile", "E5_positional_auroc"):
        print("  %-28s %10.3f %10.3f" % (k, m_raw[k], m_strat[k]))
    print("  E1 ci95  raw %s   stratified %s"
          % ([round(x, 3) for x in m_raw["E1_ci95"]],
             [round(x, 3) for x in m_strat["E1_ci95"]]))
    print("  E3 ci95  raw %s   stratified %s"
          % ([round(x, 3) for x in m_raw["E3_ci95"]],
             [round(x, 3) for x in m_strat["E3_ci95"]]))
    print("\nReference ligands under the stratified score:")
    for r in ref_strat[:12]:
        print("  %-28s arom=%d  raw %.4f  strat %.4f  %s" % (
            r["id"], n_aromatic(r["window"]), r["raw"],
            r["stratified"] if r["stratified"] is not None else float("nan"),
            r["set"]))
    print("\nTop 15 library peptides, composition-stratified:")
    for i, r in enumerate(lib_strat[:15], 1):
        print("  %2d %-14s %-26s arom=%d raw %.4f strat %.4f" % (
            i, r["peptide_id"], r["sequence"][:26], r["n_aromatic"],
            r["raw"], r["stratified"]))
    print("\nwrote %s" % out_tsv)


if __name__ == "__main__":
    main()
