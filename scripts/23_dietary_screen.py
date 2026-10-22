#!/usr/bin/env python3
"""Motif screen against the known dietary opioid peptides, with measured IC50.

This is a separate evaluation from the ChEMBL benchmark, and a better one in three
respects: the potencies are externally measured rather than curated by us, the peptides
are dietary - the actual domain of the digestion library - and the set contains both
agonists and antagonists, which the ChEMBL set does not.

The set: 20 canonical sequences from `data/dietary/dietary_peptides.tsv`, 19 with a
published IC50 spanning 0.2-200 uM (3.00 log units), 14 agonists and 5 antagonists.
CTOP is excluded as non-canonical (D-Phe, Orn, Pen).

Two pre-registered tests, both stated here before the scores were computed
----------------------------------------------------------------------------
**E8 agonist vs antagonist.** The template comes from an ACTIVE-state, agonist-bound
structure (8F7Q, mOR-beta-endorphin-Gi). Antagonists stabilise the inactive state and
bind it differently. Prediction: agonists score above antagonists.

  This test is unusually robust to the aromatic confound, because the confound points the
  OTHER WAY: the antagonists here are MORE aromatic than the agonists (mean 2.6 vs 2.0
  aromatic residues in the scored window). A score that merely counted aromatic groups
  would therefore rank antagonists HIGHER. An agonist > antagonist result cannot be
  explained by aromatic content, which is what makes it worth running.

**E9 potency ranking against measured IC50.** Spearman rho between score and pIC50 over
the peptides with a published value. Prediction: weak positive at best. The ChEMBL
equivalent (E2) gave 0.496, and this set is harder - see the caveats.

Caveats that are applied in code, not just mentioned
-----------------------------------------------------
- **Assay heterogeneity.** GPI (guinea-pig ileum) is mu-dominant; **MVD (mouse vas
  deferens) is delta-dominant**, so the two alphas1-casein exorphins were not measured at
  mu at all. They are reported but excluded from the mu-restricted analysis.
- **Terminal chemistry.** Three peptides are methyl esters in their literature form
  (casoxin-6 SRYPSY-OCH3, lactoferroxin A YLGSGY-OCH3) and one bovine lactoferrin
  fragment is not a published peptide. Esters cannot be built by this pipeline and are
  scored as free acids; all are flagged and excluded from the primary analysis.
- **IC50 is not Ki.** These are functional assays across different tissues and decades;
  the workbook itself records disputed and corrected values. Treat 3 log units of IC50 as
  an ordering, not a calibration.
- **Window truncation.** Casoxin C (YIPIQYVLSR, 10 residues) is the only peptide longer
  than the 8-residue window; everything else is scored whole.

Usage:
    python scripts/23_dietary_screen.py --runs results/dietary --out docs
"""
import argparse
import csv
import glob
import json
import random
from collections import defaultdict
from itertools import combinations
from datetime import date
from pathlib import Path

METRIC = "tversky_ref_mean"
AROMATIC = "FWYH"


def n_aromatic(w):
    return sum(w.count(c) for c in AROMATIC)


def auroc(pos, neg):
    if not pos or not neg:
        return float("nan")
    wins = sum(1.0 if p > n else (0.5 if p == n else 0.0) for p in pos for n in neg)
    return wins / (len(pos) * len(neg))


def exact_auroc_p(pos, neg):
    """One-sided exact Mann-Whitney p for AUROC(pos, neg) >= observed.

    With 3 antagonists the bootstrap CI is coarse - it can only take a few
    discrete values - so the exact permutation p is the honest statistic here.
    Enumerates every assignment of which (n_neg) of the N ranks are negatives,
    which is at most a few thousand orderings at these sizes.
    """
    npos, nneg = len(pos), len(neg)
    if not npos or not nneg:
        return float("nan")
    tot = npos * nneg
    u_obs = auroc(pos, neg) * tot
    hits = n = 0
    for c in combinations(range(1, npos + nneg + 1), nneg):
        u_pos = tot - (sum(c) - nneg * (nneg + 1) // 2)
        n += 1
        hits += u_pos >= u_obs - 1e-9
    return hits / n


def loo_antagonist_auroc(pos, neg):
    """AUROC recomputed with each negative dropped in turn.

    E8 rests on 3 antagonists. If the result is carried by one peptide, the
    range returned here will straddle 0.5 and the claim should not be made.
    """
    if len(neg) < 2:
        return []
    return sorted(auroc(pos, [x for j, x in enumerate(neg) if j != i])
                  for i in range(len(neg)))


def bootstrap_auroc(pos, neg, n=2000, seed=7):
    if not pos or not neg:
        return float("nan"), float("nan")
    rng = random.Random(seed)
    vals = sorted(auroc([rng.choice(pos) for _ in pos],
                        [rng.choice(neg) for _ in neg]) for _ in range(n))
    return vals[int(0.025 * n)], vals[int(0.975 * n) - 1]


def _ranks(xs):
    o = sorted(range(len(xs)), key=lambda i: xs[i])
    r = [0.0] * len(xs)
    i = 0
    while i < len(o):
        j = i
        while j + 1 < len(o) and xs[o[j + 1]] == xs[o[i]]:
            j += 1
        for k in range(i, j + 1):
            r[o[k]] = (i + j) / 2.0 + 1.0
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


def median(xs):
    s = sorted(xs)
    n = len(s)
    if not n:
        return float("nan")
    return s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2.0


def load_run(d):
    rows = {}
    for f in sorted(glob.glob(str(Path(d) / "block*.tsv"))):
        with open(f) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                rows[r["id"]] = r
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default=Path("results/dietary"), type=Path)
    ap.add_argument("--meta", default=Path("data/dietary/dietary_peptides.tsv"),
                    type=Path)
    ap.add_argument("--out", default=Path("docs"), type=Path)
    ap.add_argument("--jitter", default=None, type=float,
                    help="override the noise floor. By default the jitter is MEASURED "
                         "per configuration from the two seeds in that configuration's "
                         "own runs, which matters here: exclusion volumes compress the "
                         "score range (0.228 -> 0.154 median), so a jitter carried over "
                         "from a non-exclusion run is the wrong denominator. Measured, "
                         "the absolute jitter turns out NOT to compress with the scores "
                         "(0.00105 both with and without exclusion), but that is a "
                         "result, not an assumption to build in.")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    with args.meta.open() as fh:
        meta = {r["id"]: r for r in csv.DictReader(fh, delimiter="\t")}

    # configurations: <tag>_s<seed>
    runs = defaultdict(dict)
    for d in sorted(args.runs.glob("*_s*")):
        tag, _, seed = d.name.rpartition("_s")
        rows = load_run(d)
        if rows:
            runs[tag][seed] = rows
    if not runs:
        raise SystemExit(f"no dietary runs found under {args.runs}")

    report = {
        "generated": date.today().isoformat(), "metric": METRIC,
        "seed_jitter_override": args.jitter,
        "n_peptides": len(meta),
        "predictions_registered_before_scoring": {
            "E8": "agonists score above antagonists (template is an active-state, "
                  "agonist-bound structure); robust to the aromatic confound because "
                  "antagonists here are MORE aromatic (2.6 vs 2.0)",
            "E9": "weak positive correlation with pIC50 at best",
        },
        "configurations": {},
    }

    for tag in sorted(runs):
        seeds = sorted(runs[tag])
        # average the two seeds where both exist, so no single seed drives a result
        score = {}
        for pid in meta:
            vals = [float(runs[tag][s][pid][METRIC])
                    for s in seeds if pid in runs[tag][s]]
            if vals:
                score[pid] = sum(vals) / len(vals)
        if not score:
            continue

        # seed reproducibility on this set
        s1 = float("nan")
        if len(seeds) >= 2:
            shared = [p for p in meta
                      if p in runs[tag][seeds[0]] and p in runs[tag][seeds[1]]]
            if len(shared) >= 3:
                s1 = spearman([float(runs[tag][seeds[0]][p][METRIC]) for p in shared],
                              [float(runs[tag][seeds[1]][p][METRIC]) for p in shared])

        # noise floor, measured in THIS configuration rather than inherited
        jitter = args.jitter
        if jitter is None:
            jitter = float("nan")
            if len(seeds) >= 2 and len(shared) >= 3:
                jitter = median([abs(float(runs[tag][seeds[0]][p][METRIC])
                                     - float(runs[tag][seeds[1]][p][METRIC]))
                                 for p in shared])

        def subset(pred):
            return [p for p in score if pred(meta[p])]

        clean = set(subset(lambda m: m["published"] == "yes"
                           and m["terminal_supported"] == "yes"))
        mu = set(subset(lambda m: m["assay_receptor_bias"] == "mu_dominant")) & clean

        def e8(pool, label):
            ago = [score[p] for p in pool
                   if meta[p]["activity"].lower().startswith("ago")]
            ant = [score[p] for p in pool
                   if meta[p]["activity"].lower().startswith("ant")]
            if not ago or not ant:
                return None
            lo, hi = bootstrap_auroc(ago, ant)
            marg = median(ago) - median(ant)
            return {
                "pool": label, "n_agonist": len(ago), "n_antagonist": len(ant),
                "auroc": auroc(ago, ant), "ci95": [lo, hi],
                "exact_p_one_sided": exact_auroc_p(ago, ant),
                "loo_antagonist_auroc": loo_antagonist_auroc(ago, ant),
                "median_agonist": median(ago), "median_antagonist": median(ant),
                "margin": marg, "seed_jitter_measured": jitter,
                "margin_in_jitter": marg / jitter,
                "noise_floor_limited": bool(abs(marg) < 2 * jitter),
                "mean_aromatic_agonist": sum(
                    n_aromatic(meta[p]["sequence"][:8]) for p in pool
                    if meta[p]["activity"].lower().startswith("ago")) / max(1, len(ago)),
                "mean_aromatic_antagonist": sum(
                    n_aromatic(meta[p]["sequence"][:8]) for p in pool
                    if meta[p]["activity"].lower().startswith("ant")) / max(1, len(ant)),
            }

        def e9(pool, label, agonists_only=False):
            pts = [(score[p], float(meta[p]["pIC50"])) for p in pool
                   if meta[p]["pIC50"]
                   and (not agonists_only
                        or meta[p]["activity"].lower().startswith("ago"))]
            if len(pts) < 4:
                return None
            return {"pool": label, "n": len(pts),
                    "spearman_vs_pIC50": spearman([a for a, _ in pts],
                                                  [b for _, b in pts])}

        # confound check on this set
        arom = [(score[p], n_aromatic(meta[p]["sequence"][:8])) for p in clean]
        report["configurations"][tag] = {
            "seeds": seeds, "n_scored": len(score),
            "S1_spearman_this_set": s1,
            "E8_all_published": e8(clean, "published, supported termini"),
            "E8_mu_assays_only": e8(mu, "mu-dominant assays only"),
            "E9_all_published": e9(clean, "published, supported termini"),
            "E9_mu_assays_only": e9(mu, "mu-dominant assays only"),
            "E9_mu_agonists_only": e9(mu, "mu-dominant agonists only",
                                      agonists_only=True),
            "E6_rho_vs_aromatic": spearman([a for a, _ in arom],
                                           [b for _, b in arom]),
            "ranked": sorted(
                ({"id": p, "name": meta[p]["name"], "sequence": meta[p]["sequence"],
                  "activity": meta[p]["activity"], "ic50_uM": meta[p]["ic50_uM"],
                  "pIC50": meta[p]["pIC50"], "score": score[p],
                  "n_aromatic": n_aromatic(meta[p]["sequence"][:8]),
                  "assay_bias": meta[p]["assay_receptor_bias"],
                  "flagged": meta[p]["published"] != "yes"
                             or meta[p]["terminal_supported"] != "yes"}
                 for p in score), key=lambda r: -r["score"]),
        }

    (args.out / "dietary_screen.json").write_text(json.dumps(report, indent=2) + "\n")

    for tag, c in report["configurations"].items():
        print("=" * 72)
        print("template %s   seeds %s   S1 on this set %.3f"
              % (tag, ",".join(c["seeds"]), c["S1_spearman_this_set"]))
        print("  E6 rho vs aromatic count on this set: %+.3f" % c["E6_rho_vs_aromatic"])
        for key in ("E8_all_published", "E8_mu_assays_only"):
            e = c[key]
            if not e:
                continue
            print("  %-20s AUROC %.3f [%.2f, %.2f]  n=%d/%d  margin %+.4f (%+.1fx "
                  "jitter, measured %.5f)%s" % (key, e["auroc"], e["ci95"][0], e["ci95"][1],
                                 e["n_agonist"], e["n_antagonist"], e["margin"],
                                 e["margin_in_jitter"], e["seed_jitter_measured"],
                                 "  NOISE-FLOOR" if e["noise_floor_limited"] else ""))
            print("  %-20s exact p %.4f (one-sided)   leave-one-antagonist-out "
                  "AUROC %s" % ("", e["exact_p_one_sided"],
                                " ".join("%.3f" % v
                                         for v in e["loo_antagonist_auroc"])))
            print("  %-20s aromatic: agonists %.2f vs antagonists %.2f"
                  % ("", e["mean_aromatic_agonist"], e["mean_aromatic_antagonist"]))
        for key in ("E9_all_published", "E9_mu_assays_only", "E9_mu_agonists_only"):
            e = c[key]
            if e:
                print("  %-20s rho %+.3f  (n=%d)" % (key, e["spearman_vs_pIC50"],
                                                     e["n"]))
        print("\n  rank  name                           seq        act   IC50    arom"
              "  score")
        for i, r in enumerate(c["ranked"], 1):
            print("  %4d  %-30s %-10s %-5s %-7s %d     %.4f%s" % (
                i, r["name"][:30], r["sequence"][:10], r["activity"][:3],
                r["ic50_uM"] or "-", r["n_aromatic"], r["score"],
                "  (flagged)" if r["flagged"] else ""))
    print("\nwrote %s" % (args.out / "dietary_screen.json"))


if __name__ == "__main__":
    main()
