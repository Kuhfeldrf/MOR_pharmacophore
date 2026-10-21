#!/usr/bin/env python3
"""Every AUROC in this study, restated as an effect size against the noise floor.

Why this exists
---------------
An AUROC says how consistently two groups are ordered. It says nothing about how far
apart they are. This project has now been misled twice by that gap:

  - The first study reported a tier structure that was 73% decile churn between seeds.
  - This study reported E1 = 1.000 at 400 conformers, which fell to 0.803 when the
    ensemble was doubled and the ordering stopped being an artefact of under-sampling.

Both would have been caught immediately by reporting the separation in units of the
measured seed jitter. So that is now reported for every pre-registered contrast:

    margin = median(positive scores) - median(negative scores)
    margin_in_jitter = margin / S1_median_abs_delta

A contrast whose margin is below about 2x the jitter is NOISE-FLOOR LIMITED: its AUROC
may be well above chance and still describe an ordering of differences the method cannot
reliably reproduce. That is not the same as being wrong, and it is not the same as being
usable. It is reported as its own status so it cannot be read as either.

This is a reporting change, not a new criterion for choosing anything. It is applied to
numbers that were already computed, and it changes no selection.

Usage:
    python scripts/18_effect_sizes.py --bench-run results/stageB/nc800_w8_s61453 \\
        --epimers results/epimers --out docs
"""
import argparse
import csv
import glob
import json
from collections import defaultdict
from datetime import date
from pathlib import Path

METRIC = "tversky_ref_mean"
NOISE_FLOOR_MULTIPLE = 2.0


def median(xs):
    if not xs:
        return float("nan")
    s = sorted(xs)
    n = len(s)
    return s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2.0


def auroc(pos, neg):
    if not pos or not neg:
        return float("nan")
    wins = sum(1.0 if p > n else (0.5 if p == n else 0.0) for p in pos for n in neg)
    return wins / (len(pos) * len(neg))


def load(run, key=METRIC):
    out = {}
    for f in sorted(glob.glob(str(Path(run) / "block*.tsv"))):
        with open(f) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                out[r["id"]] = r
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bench-run", default=Path("results/stageB/nc800_w8_s61453"),
                    type=Path)
    ap.add_argument("--epimers", default=Path("results/epimers"), type=Path)
    ap.add_argument("--benchmark", default=Path("data/benchmark/benchmark.tsv"),
                    type=Path)
    ap.add_argument("--sweep", default=Path("docs/sweep_results.json"), type=Path)
    ap.add_argument("--out", default=Path("docs"), type=Path)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    # the measured noise floor for the gated configuration
    jitter = None
    for r in json.loads(args.sweep.read_text())["results"]:
        if r["window"] == 8 and r["n_conformers"] == 800:
            jitter = r["S1_median_abs_delta"]
    if jitter is None:
        raise SystemExit("could not read S1_median_abs_delta for window 8 / 800")

    rows = load(args.bench_run)
    with args.benchmark.open() as fh:
        bench = {r["id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    sets, fam = defaultdict(list), defaultdict(list)
    for bid, r in bench.items():
        if bid not in rows:
            continue
        v = float(rows[bid][METRIC])
        sets[r["set"]].append(v)
        if r["set"] == "active" and r["family"]:
            fam[r["family"]].append(v)

    contrasts = []

    def add(name, pos, neg, desc):
        m = median(pos) - median(neg)
        contrasts.append({
            "contrast": name, "description": desc,
            "auroc": auroc(pos, neg),
            "n_pos": len(pos), "n_neg": len(neg),
            "median_pos": median(pos), "median_neg": median(neg),
            "margin": m, "margin_in_jitter": m / jitter if jitter else float("nan"),
            "noise_floor_limited": bool(abs(m) < NOISE_FLOOR_MULTIPLE * jitter),
        })

    add("S2", sets["active"], sets["scramble"], "actives vs composition-matched scrambles")
    add("E1", fam.get("YGGFL", []), fam.get("YGPFL", []),
        "YGGFL family vs YGPFL family (Gly3->Pro, ~3 log units)")
    add("E3", sets["active"], sets["tyr_hard"],
        "actives vs Tyr-starting length-matched library decoys")
    add("E5", sets["positional_start"], sets["positional_span"],
        "starts-with-YPFP vs spans-YPFP")

    # E7, from the epimer run: paired, so the margin is the median paired difference
    epi = load(args.epimers)
    pair_file = Path("data/epimers/epimer_pairs.tsv")
    if epi and pair_file.exists():
        with pair_file.open() as fh:
            meta = {r["id"]: r for r in csv.DictReader(fh, delimiter="\t")}
        pairs = defaultdict(dict)
        for pid, m in meta.items():
            if pid in epi:
                pairs[m["pair"]][m["form"]] = float(epi[pid][METRIC])
        mu = ("dermorphin", "dala2_met_enkephalin", "dala2_leu_enkephalin")
        deltas = [v["D"] - v["L"] for k, v in pairs.items()
                  if k in mu and "D" in v and "L" in v]
        if deltas:
            m = median(deltas)
            contrasts.append({
                "contrast": "E7", "description":
                    "mu-relevant L/D epimer pairs, paired difference (exactly "
                    "composition-matched)",
                "auroc": None, "n_pos": len(deltas), "n_neg": len(deltas),
                "median_pos": None, "median_neg": None,
                "margin": m, "margin_in_jitter": m / jitter,
                "noise_floor_limited": bool(abs(m) < NOISE_FLOOR_MULTIPLE * jitter),
                "paired_deltas": deltas,
                "n_correct_direction": sum(1 for d in deltas if d > 0),
            })

    payload = {
        "generated": date.today().isoformat(),
        "metric": METRIC,
        "configuration": {"window": 8, "n_conformers": 800},
        "seed_jitter_median_abs_delta": jitter,
        "noise_floor_multiple": NOISE_FLOOR_MULTIPLE,
        "note": ("Reporting change only; applied to already-computed numbers and used "
                 "to choose nothing. A contrast below the noise-floor multiple has an "
                 "ordering the method cannot reliably reproduce, whatever its AUROC."),
        "contrasts": contrasts,
    }
    (args.out / "effect_sizes.json").write_text(json.dumps(payload, indent=2) + "\n")

    print("Effect sizes against the measured noise floor")
    print("window 8 / 800 conformers, %s" % METRIC)
    print("seed jitter (S1 median |delta|) = %.4f; noise-floor threshold = %.1fx = "
          "%.4f\n" % (jitter, NOISE_FLOOR_MULTIPLE, NOISE_FLOOR_MULTIPLE * jitter))
    print("  contrast  AUROC   median+   median-    margin   in jitter  status")
    for c in contrasts:
        a = "  n/a " if c["auroc"] is None else "%6.3f" % c["auroc"]
        mp = "   n/a " if c["median_pos"] is None else "%.4f" % c["median_pos"]
        mn = "   n/a " if c["median_neg"] is None else "%.4f" % c["median_neg"]
        print("  %-8s %s  %s   %s  %+.4f   %+6.2fx   %s" % (
            c["contrast"], a, mp, mn, c["margin"], c["margin_in_jitter"],
            "NOISE-FLOOR LIMITED" if c["noise_floor_limited"] else "resolved"))
    print()
    for c in contrasts:
        print("  %-4s %s" % (c["contrast"], c["description"]))
    print("\nwrote %s" % (args.out / "effect_sizes.json"))


if __name__ == "__main__":
    main()
