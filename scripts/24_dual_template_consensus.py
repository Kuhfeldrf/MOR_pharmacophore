#!/usr/bin/env python3
"""Dual-template analysis: consensus, paired template comparison, and the re-pointed E4'.

Design and pre-registered predictions: `docs/dual_template_design.md`. This script
computes the tests in its section 6 and nothing else; it does not choose anything, and
every threshold it checks against was written down before the scores existed.

Division of labour with `23_dietary_screen.py`
----------------------------------------------
That script already discovers every `<tag>_s<seed>` directory under `results/dietary` and
reports E8 (agonist vs antagonist AUROC) and E9 (potency ranking) per configuration, so
the new `b8excl` and `b4excl` runs get their single-template numbers from it for free, on
exactly the same code path as the existing `p8excl` runs. Running it twice would be two
implementations of one statistic, which is how two numbers that should agree stop
agreeing. Its helpers are imported here rather than reimplemented.

What is genuinely new, and therefore lives here:

  T2  the consensus statistic over two templates (rank fusion, and max-normalised)
  T3  the paired comparison at identical window, including rho between the two
      templates' score vectors - the test of whether template B adds information or
      merely re-measures template A
  T5  E4', the held-out DAMGO percentile against the library distribution
  T6  rho(score, aromatic count) compared between templates

Usage:
    python scripts/24_dual_template_consensus.py \\
        --runs results/dietary --template-a p8excl --template-b b8excl \\
        --holdout results/holdout --library results/library_b8excl --out docs
"""
import argparse
import csv
import glob
import importlib.util
import json
from datetime import date
from pathlib import Path

_HERE = Path(__file__).resolve().parent


def _dietary_helpers():
    """Import `23_dietary_screen.py`, whose name is not a legal module identifier."""
    spec = importlib.util.spec_from_file_location(
        "dietary_screen", _HERE / "23_dietary_screen.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


D = _dietary_helpers()
METRIC = D.METRIC

# Pre-registered in docs/dual_template_design.md 6, before any score was computed.
PREREGISTERED = {
    "T1": "B8 agonist-vs-antagonist AUROC >= 0.75 and point estimate >= A8",
    "T2": "consensus AUROC >= max(A8, B8)",
    "T3": "B8 >= A8 at identical window, and rho(A8, B8) < 0.9",
    "T4": "B8 Spearman vs pIC50 in (0, 0.6)",
    "T5": "DAMGO at or above the 90th percentile of the library for >= 1 template",
    "T6": "rho(score, aromatic count) for B8 < 0.909",
}


def classify(activity):
    a = (activity or "").strip().lower()
    if "antagonist" in a:
        return "antagonist"
    if "agonist" in a:
        return "agonist"
    return "other"


def rank_fuse(score_by_template, ids):
    """Mean rank across templates. Higher is better.

    Ranks rather than raw scores because the two templates' score distributions are not
    on a common scale: template A carries 27 features from 8 residues, template B carries
    the features of a complete 4-residue ligand, and Tversky reference similarity is
    normalised by the reference, so the same peptide-pocket fit gives systematically
    different numbers. Rank fusion needs no calibration between them; a raw average would
    silently weight whichever template happens to have the wider spread.
    """
    per_template_ranks = []
    for scores in score_by_template:
        ordered = sorted(ids, key=lambda i: scores[i])
        rank = {i: n for n, i in enumerate(ordered, start=1)}
        per_template_ranks.append(rank)
    return {i: sum(r[i] for r in per_template_ranks) / len(per_template_ranks)
            for i in ids}


def max_normalised(score_by_template, ids):
    """Per-peptide maximum of min-max scaled scores.

    The biologically motivated reading: a dietary peptide is interesting if it mimics
    EITHER endogenous binding mode, so it should not be penalised for failing to match
    the mode it was never going to match. Reported as secondary because max-of-two is
    upward-biased by construction - two chances to score well - and the antagonists get
    that same bias, so it does not cancel out.
    """
    scaled = []
    for scores in score_by_template:
        vals = [scores[i] for i in ids]
        lo, hi = min(vals), max(vals)
        rng = (hi - lo) or 1.0
        scaled.append({i: (scores[i] - lo) / rng for i in ids})
    return {i: max(s[i] for s in scaled) for i in ids}


def pools(ids, meta):
    """The two E8 pools, defined exactly as `23_dietary_screen.py` defines them.

    This is not a detail. E8 does not score all 20 dietary peptides: the published 0.917
    is the mu-dominant-assay pool, 12 agonists against 3 antagonists. Computing an AUROC
    over all 20 instead gives a different number that looks like it should be comparable
    and is not - the two alphas1-casein exorphins were measured on MVD and never at mu at
    all, and three peptides have unsupported terminal chemistry. Duplicating the pool
    definition here rather than importing it would be the same hazard one level down, so
    the predicates are copied verbatim from that script and the pool labels are kept
    identical.
    """
    clean = {i for i in ids
             if meta[i]["published"] == "yes" and meta[i]["terminal_supported"] == "yes"}
    mu = {i for i in clean if meta[i]["assay_receptor_bias"] == "mu_dominant"}
    return [("mu-dominant assays only", sorted(mu)),
            ("published, supported termini", sorted(clean))]


def auroc_block(scores, ids, meta):
    ago = [scores[i] for i in ids if classify(meta[i]["activity"]) == "agonist"]
    ant = [scores[i] for i in ids if classify(meta[i]["activity"]) == "antagonist"]
    if not ago or not ant:
        return None
    lo, hi = D.bootstrap_auroc(ago, ant)
    return {
        "n_agonist": len(ago), "n_antagonist": len(ant),
        "auroc": D.auroc(ago, ant),
        "ci95": [lo, hi],
        "exact_p_one_sided": D.exact_auroc_p(ago, ant),
        # Full leave-one-antagonist-out range, not a summary: with only 3 antagonists the
        # question is whether the result survives dropping any single one, and a range
        # that straddles 0.5 means it does not (D.loo_antagonist_auroc docstring).
        "loo_antagonist_auroc": D.loo_antagonist_auroc(ago, ant),
        "median_agonist": D.median(ago), "median_antagonist": D.median(ant),
    }


def by_pool(scores, ids, meta):
    """`auroc_block` for each E8 pool, keyed by the pool's label."""
    out = {}
    for label, pool in pools(ids, meta):
        pool = [i for i in pool if i in scores]
        blk = auroc_block(scores, pool, meta) if pool else None
        if blk:
            out[label] = blk
    return out


def percentile_of(value, population):
    """Fraction of the population at or below `value`, as a percentile."""
    if not population:
        return float("nan")
    return 100.0 * sum(1 for p in population if p <= value) / len(population)


def load_library_scores(d):
    """Every library peptide's metric value under one template."""
    vals = []
    for f in sorted(glob.glob(str(Path(d) / "block*.tsv"))):
        with open(f) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("set") != "library":
                    continue          # reference ligands must not be in their own baseline
                try:
                    vals.append(float(r[METRIC]))
                except (KeyError, ValueError):
                    continue
    return vals


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default=Path("results/dietary"), type=Path)
    ap.add_argument("--meta", default=Path("data/dietary/dietary_peptides.tsv"), type=Path)
    ap.add_argument("--template-a", default="p8excl",
                    help="tag of the 8F7Q/beta-endorphin run (window 8, exclusion)")
    ap.add_argument("--template-b", default="b8excl",
                    help="tag of the 8F7R/endomorphin-1 run (window 8, exclusion)")
    ap.add_argument("--holdout", default=Path("results/holdout"), type=Path)
    ap.add_argument("--library", default=None, type=Path,
                    help="library run under template B; needed for T5 only. Absent "
                         "means T5 is reported as not yet run, not as a failure.")
    ap.add_argument("--out", default=Path("docs"), type=Path)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    with args.meta.open() as fh:
        meta = {r["id"]: r for r in csv.DictReader(fh, delimiter="\t")}

    # Discover <tag>_s<seed> exactly as 23_dietary_screen.py does.
    runs = {}
    for d in sorted(args.runs.glob("*_s*")):
        tag, _, seed = d.name.rpartition("_s")
        rows = D.load_run(d)
        if rows:
            runs.setdefault(tag, {})[seed] = rows

    for tag in (args.template_a, args.template_b):
        if tag not in runs:
            raise SystemExit(
                f"no runs for template tag {tag!r} under {args.runs}. "
                f"Found: {sorted(runs)}. Run jobs/dual_template.sbatch first.")

    seeds = sorted(set(runs[args.template_a]) & set(runs[args.template_b]))
    if not seeds:
        raise SystemExit(
            f"templates {args.template_a!r} and {args.template_b!r} share no seed; "
            "a paired comparison across different seeds would confound template with "
            "conformer ensemble, which is the one thing T3 exists to separate.")

    report = {
        "generated": date.today().isoformat(),
        "metric": METRIC,
        "design": "docs/dual_template_design.md",
        "preregistered": PREREGISTERED,
        "template_a": {"tag": args.template_a, "source": "8F7Q", "ligand": "beta-endorphin(1-31)"},
        "template_b": {"tag": args.template_b, "source": "8F7R", "ligand": "endomorphin-1"},
        "seeds": seeds,
        "per_seed": {},
    }

    def analyse(a, b, ids, windows):
        """One comparison block from two aligned score maps."""
        fused = rank_fuse([a, b], ids)
        maxnorm = max_normalised([a, b], ids)
        arom = {i: D.n_aromatic(windows.get(i) or meta[i]["sequence"][:8]) for i in ids}
        # Potency subset: mu-dominant assays with supported terminal chemistry only. The
        # two alphas1-casein exorphins were measured on MVD (delta-dominant) and so were
        # never measured at mu; esters are scored as free acids and are excluded too.
        pot = [(i, float(meta[i]["pIC50"])) for i in ids
               if meta[i].get("pIC50")
               and meta[i].get("assay_receptor_bias") == "mu_dominant"
               and meta[i].get("terminal_supported") == "yes"]
        return {
            "n_scored": len(ids),
            "T1_template_a": by_pool(a, ids, meta),
            "T1_template_b": by_pool(b, ids, meta),
            "T2_consensus_rank_fusion": by_pool(fused, ids, meta),
            "T2_consensus_max_normalised": by_pool(maxnorm, ids, meta),
            "T3_rho_between_templates": D.spearman([a[i] for i in ids], [b[i] for i in ids]),
            "T4_spearman_vs_pIC50": {
                "template_a": D.spearman([a[i] for i, _ in pot], [p for _, p in pot]),
                "template_b": D.spearman([b[i] for i, _ in pot], [p for _, p in pot]),
                "n": len(pot),
            },
            "T6_rho_vs_aromatic": {
                "template_a": D.spearman([a[i] for i in ids], [arom[i] for i in ids]),
                "template_b": D.spearman([b[i] for i in ids], [arom[i] for i in ids]),
            },
        }

    def scores_for(tag, seed_list):
        """Mean metric per peptide across `seed_list`, plus the scored window."""
        rows_by_seed = [runs[tag][s] for s in seed_list]
        ids = set.intersection(*(set(r) for r in rows_by_seed)) & set(meta)
        vals, windows = {}, {}
        for i in ids:
            vals[i] = sum(float(r[i][METRIC]) for r in rows_by_seed) / len(rows_by_seed)
            windows[i] = rows_by_seed[0][i].get("window")
        return vals, windows

    # Primary: seed-averaged, matching 23_dietary_screen.py so T1 is directly comparable
    # with the published 0.917 rather than merely resembling it.
    a_avg, win = scores_for(args.template_a, seeds)
    b_avg, _ = scores_for(args.template_b, seeds)
    shared = sorted(set(a_avg) & set(b_avg))
    report["seed_averaged"] = analyse(a_avg, b_avg, shared, win)

    # Secondary: each seed alone, to check nothing reverses.
    for seed in seeds:
        a_one, win_one = scores_for(args.template_a, [seed])
        b_one, _ = scores_for(args.template_b, [seed])
        ids = sorted(set(a_one) & set(b_one))
        report["per_seed"][seed] = analyse(a_one, b_one, ids, win_one)

    # --- T5: E4' held-out retrieval -------------------------------------------------
    holdout = {}
    for d in sorted(Path(args.holdout).glob("*_s*")) if Path(args.holdout).exists() else []:
        tag, _, seed = d.name.rpartition("_s")
        rows = D.load_run(d)
        if rows:
            holdout.setdefault(tag, {})[seed] = rows
    library = load_library_scores(args.library) if args.library and Path(args.library).exists() else []
    report["T5_holdout"] = {
        "library_run": str(args.library) if args.library else None,
        "n_library_scored": len(library),
        "status": "computed" if library else
                  "library distribution not present - run jobs/dual_template_library.sbatch",
        "ligands": {},
    }
    for tag, by_seed in sorted(holdout.items()):
        for seed, rows in sorted(by_seed.items()):
            for lid, row in sorted(rows.items()):
                try:
                    v = float(row[METRIC])
                except (KeyError, ValueError):
                    continue
                report["T5_holdout"]["ligands"].setdefault(lid, {})[f"{tag}_s{seed}"] = {
                    "score": v,
                    "percentile_vs_library": percentile_of(v, library) if library else None,
                }

    out_json = args.out / "dual_template_result.json"
    out_json.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")

    # --- console summary -------------------------------------------------------------
    print(f"metric {METRIC}  seeds {','.join(seeds)}")
    print(f"template A {args.template_a} (8F7Q, beta-endorphin)")
    print(f"template B {args.template_b} (8F7R, endomorphin-1)\n")
    def show(label, blk):
        print(f"--- {label}  n={blk['n_scored']} ---")
        for key in ("T1_template_a", "T1_template_b",
                    "T2_consensus_rank_fusion", "T2_consensus_max_normalised"):
            for pool, e in sorted(blk[key].items()):
                loo = e["loo_antagonist_auroc"]
                loo_s = f"{min(loo):.3f}-{max(loo):.3f}" if loo else "n/a"
                straddles = loo and min(loo) < 0.5 < max(loo)
                print(f"  {key:30} [{pool:28}] n={e['n_agonist']}/{e['n_antagonist']}"
                      f" AUROC {e['auroc']:.3f} p={e['exact_p_one_sided']:.4f}"
                      f" LOO {loo_s}{'  <-- LOO straddles 0.5, not claimable' if straddles else ''}")
        rho = blk["T3_rho_between_templates"]
        print(f"  T3 rho(A,B)                    {rho:+.3f}"
              f"   {'REDUNDANT (>=0.9): second template adds nothing' if abs(rho) >= 0.9 else 'templates differ'}")
        t4 = blk["T4_spearman_vs_pIC50"]
        print(f"  T4 rho vs pIC50 (n={t4['n']})          A {t4['template_a']:+.3f}   B {t4['template_b']:+.3f}")
        t6 = blk["T6_rho_vs_aromatic"]
        print(f"  T6 rho vs aromatic count       A {t6['template_a']:+.3f}   B {t6['template_b']:+.3f}"
              f"   (prior truncated-template value +0.909)")
        print()

    show("SEED-AVERAGED (primary)", report["seed_averaged"])
    for seed, blk in report["per_seed"].items():
        show(f"seed {seed} alone", blk)
    print(f"T5 {report['T5_holdout']['status']}")
    for lid, per in sorted(report["T5_holdout"]["ligands"].items()):
        for run, v in sorted(per.items()):
            pc = v["percentile_vs_library"]
            print(f"  {lid:28} {run:18} score {v['score']:.4f}"
                  + (f"  {pc:.1f}th pct" if pc is not None else "  (no library baseline)"))
    print(f"\nwrote {out_json}")


if __name__ == "__main__":
    main()
