#!/usr/bin/env python3
"""Stage 6: evaluate the pre-registered acceptance criteria (PLAN.md 3).

C1  template recovery          the template's own ligand scores in the top decile
C2  independent-ligand recovery endomorphin, from a structure not used to build it
C3  added value over baseline  separation WITHIN the Tyr-starting subset, where the
                               "starts with Tyr" baseline is constant
C4  positional rule            peptides starting with YPFP vs those merely spanning it
C5  conformer stability        tier churn between two seeds (needs a second run)

C3 is the criterion that matters. C1, C2 and C4 can all pass while the screen still
reduces to the free sequence rule, which is the failure mode PLAN.md 2 names.

No criterion is adjusted here to fit a result. Usage:
    python scripts/05_validate.py --screen results/04_screen --out docs
"""
import argparse
import csv
import json
import random
from datetime import date
from pathlib import Path

METRIC = "tversky_ref"
TOP_DECILE = 0.90


def load_screen(screen_dir):
    rows = []
    for tsv in sorted(Path(screen_dir).glob("block*.tsv")):
        with tsv.open() as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                r[METRIC] = float(r[METRIC])
                r["tversky_ref_spread"] = float(r["tversky_ref_spread"])
                rows.append(r)
    return rows


def quantile_of(value, population):
    """Fraction of the population at or below `value`."""
    if not population:
        return float("nan")
    return sum(1 for p in population if p <= value) / len(population)


def auroc(positives, negatives):
    """Rank-based AUROC with explicit tie handling (ties contribute 0.5)."""
    if not positives or not negatives:
        return float("nan")
    wins = 0.0
    for p in positives:
        for n in negatives:
            if p > n:
                wins += 1.0
            elif p == n:
                wins += 0.5
    return wins / (len(positives) * len(negatives))


def bootstrap_auroc(positives, negatives, n=2000, seed=11):
    rng = random.Random(seed)
    vals = []
    for _ in range(n):
        p = [rng.choice(positives) for _ in positives]
        q = [rng.choice(negatives) for _ in negatives]
        vals.append(auroc(p, q))
    vals.sort()
    lo = vals[int(0.025 * len(vals))]
    hi = vals[int(0.975 * len(vals)) - 1]
    return lo, hi


def median(xs):
    if not xs:
        return float("nan")
    s = sorted(xs)
    m = len(s) // 2
    return s[m] if len(s) % 2 else (s[m - 1] + s[m]) / 2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--screen", default=Path("results/04_screen"), type=Path)
    ap.add_argument("--controls", default=Path("data/controls/control_sets.json"),
                    type=Path)
    ap.add_argument("--screen2", default=None,
                    help="a second-seed screen directory, for the C5 stability check")
    ap.add_argument("--out", default=Path("docs"), type=Path)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    rows = load_screen(args.screen)
    if not rows:
        raise SystemExit(f"no block*.tsv under {args.screen}")
    controls = json.loads(args.controls.read_text())

    by_id = {r["id"]: r for r in rows}
    by_seq = {}
    for r in rows:
        by_seq.setdefault(r["sequence"], r)

    library = [r for r in rows if r["set"] == "library"]
    library_scores = [r[METRIC] for r in library]
    results = {"generated": date.today().isoformat(), "metric": METRIC,
               "n_scored": len(rows), "n_library": len(library)}

    def score_of(sequence, ref_name=None):
        if ref_name and f"ref__{ref_name}" in by_id:
            return by_id[f"ref__{ref_name}"][METRIC]
        r = by_seq.get(sequence)
        return r[METRIC] if r else None

    # ---- C1 template recovery -------------------------------------------------
    tmpl = controls["template_ligand"][0]
    s = score_of(tmpl["sequence"], tmpl["name"])
    q = quantile_of(s, library_scores) if s is not None else float("nan")
    results["C1"] = {"ligand": tmpl["name"], "sequence": tmpl["sequence"],
                     "score": s, "library_quantile": q,
                     "pass": bool(s is not None and q >= TOP_DECILE)}

    # ---- C2 independent-ligand recovery ---------------------------------------
    endo = next((p for p in controls["positive_endogenous"]
                 if p["name"] == "endomorphin_1"), None)
    s2 = score_of(endo["sequence"], endo["name"]) if endo else None
    q2 = quantile_of(s2, library_scores) if s2 is not None else float("nan")
    results["C2"] = {"ligand": endo["name"] if endo else None,
                     "score": s2, "library_quantile": q2,
                     "pass": bool(s2 is not None and q2 >= TOP_DECILE)}

    # ---- C3 added value over the sequence baseline ----------------------------
    # Positives: documented actives. Negatives: Tyr-starting library peptides that are
    # not known actives. Restricting to Tyr-starting holds the baseline constant.
    active_seqs = {p["sequence"] for p in controls["positive_endogenous"]}
    active_seqs |= {p["sequence"] for p in controls["positive_food"]}

    pos_tyr, neg_tyr = [], []
    for p in controls["positive_endogenous"] + controls["positive_food"]:
        if not p["sequence"].startswith("Y"):
            continue
        v = score_of(p["sequence"], p["name"])
        if v is not None:
            pos_tyr.append(v)
    for r in library:
        if r["sequence"].startswith("Y") and r["sequence"] not in active_seqs:
            neg_tyr.append(r[METRIC])

    a3 = auroc(pos_tyr, neg_tyr)
    lo3, hi3 = bootstrap_auroc(pos_tyr, neg_tyr) if pos_tyr and neg_tyr else (float("nan"),) * 2
    # Effective n: peptides sharing an N-terminal window are one entity to this screen.
    windows_neg = {r["sequence"][:5] for r in library
                   if r["sequence"].startswith("Y") and r["sequence"] not in active_seqs}
    results["C3"] = {
        "n_positives": len(pos_tyr), "n_negatives": len(neg_tyr),
        "effective_n_negatives_unique_windows": len(windows_neg),
        "auroc": a3, "ci95": [lo3, hi3],
        "pass": bool(a3 == a3 and a3 >= 0.70 and lo3 > 0.5),
        "note": ("Negatives share N-terminal windows, so ties are pervasive and the "
                 "effective n is the unique-window count, not the peptide count."),
    }

    # ---- C4 positional rule ----------------------------------------------------
    pc = controls["positional_control"]
    starts = [by_id[e["peptide_id"]][METRIC] for e in pc["starts_with_motif"]
              if e["peptide_id"] in by_id]
    spans = [by_id[e["peptide_id"]][METRIC] for e in pc["spans_motif_blocked_nterm"]
             if e["peptide_id"] in by_id]
    a4 = auroc(starts, spans)
    results["C4"] = {
        "motif": pc["motif"], "n_starts": len(starts), "n_spans": len(spans),
        "median_starts": median(starts), "median_spans": median(spans),
        "median_difference": (median(starts) - median(spans))
        if starts and spans else float("nan"),
        "auroc": a4,
        "pass": bool(a4 == a4 and a4 > 0.5),
    }

    # ---- C5 conformer stability ------------------------------------------------
    spreads = [r["tversky_ref_spread"] for r in rows]
    c5 = {
        "median_ensemble_spread": median(spreads),
        "max_ensemble_spread": max(spreads) if spreads else float("nan"),
    }
    if args.screen2 and Path(args.screen2).exists():
        rows2 = load_screen(args.screen2)
        other = {r["id"]: r[METRIC] for r in rows2}
        shared = [r for r in rows if r["id"] in other]
        # Tiering is defined by decile of the library distribution, so "tier churn" is
        # the fraction of peptides that move decile between seeds.
        lib2 = sorted(other[r["id"]] for r in rows2 if r["set"] == "library")

        def decile(value, population):
            if not population:
                return -1
            return min(9, int(quantile_of(value, population) * 10))

        churn = sum(
            1 for r in shared
            if decile(r[METRIC], library_scores) != decile(other[r["id"]], lib2))
        frac = churn / len(shared) if shared else float("nan")
        deltas = [abs(r[METRIC] - other[r["id"]]) for r in shared]
        c5.update({
            "seed2": str(args.screen2), "n_compared": len(shared),
            "decile_churn_fraction": frac,
            "median_abs_score_delta": median(deltas),
            "max_abs_score_delta": max(deltas) if deltas else float("nan"),
            "pass": bool(frac == frac and frac <= 0.10),
        })
    else:
        c5["status"] = "no second-seed run supplied; ensemble spread only"
        c5["pass"] = None
    results["C5"] = c5

    results["library_score_distribution"] = {
        "min": min(library_scores), "median": median(library_scores),
        "max": max(library_scores),
        "p90": sorted(library_scores)[int(0.90 * len(library_scores))],
    }

    (args.out / "screen_validation.json").write_text(json.dumps(results, indent=2) + "\n")

    for key in ("C1", "C2", "C3", "C4", "C5"):
        print(f"{key}: {json.dumps(results[key])}")
    print(f"\nwrote {args.out / 'screen_validation.json'}")


if __name__ == "__main__":
    main()
