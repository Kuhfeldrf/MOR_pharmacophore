#!/usr/bin/env python3
"""Milk-atlas screen analysis, with predictions registered before the scores exist.

The atlas set is not the dietary set and cannot reuse `23_dietary_screen.py`: there are no
agonist/antagonist labels and no measured IC50 here, because four of the five novel
candidates have never been assayed. What this set does have is **five released sites whose
peptides are published opioids** (beta-casomorphin-7, human beta-casomorphin-5,
neocasomorphin-6, the alphas1-casein exorphin, and the conserved alpha-lactorphin site) and
**five with no published activity**. So the available test is whether the screen puts the
known ones first.

Predictions, written down before any atlas score was computed (docs/milk_atlas.md 5):

  M1  PRIMARY. Known-opioid sites rank above novel sites. AUROC(known, novel) >= 0.70.
      This is the positive-control test, and it is the one that licenses reading
      anything into the novel scores at all. n = 5 vs 5, so it is underpowered by
      construction and reported with an exact permutation p and a
      leave-one-known-out range.

  M2  The two aromatic-rich novel candidates (YPYAVWYYPQ, YPTYHTHGRYVPP; 4 aromatics
      each) lose MORE score when exclusion volumes are applied than the known controls
      do. Exclusion volumes encode the chamber's steric bound, and the failure mode they
      exist to catch is a peptide that scores well by over-filling the pocket with
      aromatic groups. If these two behave like the known opioids under exclusion, that
      is evidence they are real candidates rather than confound artefacts.

  M3  rho(score, aromatic count) is weaker WITH exclusion volumes than without, on this
      set. The prior value on the benchmark set was +0.909 without.

  M4  Seed stability: Spearman rho between the two seeds >= 0.95, the S1 gate. A
      configuration failing this cannot support a shortlist.

Nothing here measures binding. A peptide at the top of this ranking is a synthesis
candidate, not a hit.

Usage:
    python scripts/27_milk_atlas_analysis.py --runs results/milk_atlas \\
        --targets data/milk_atlas/milk_atlas_targets.tsv --out docs
"""
import argparse
import csv
import importlib.util
import json
from collections import defaultdict
from datetime import date
from pathlib import Path

_HERE = Path(__file__).resolve().parent


def _dietary_helpers():
    spec = importlib.util.spec_from_file_location(
        "dietary_screen", _HERE / "23_dietary_screen.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


D = _dietary_helpers()
METRIC = D.METRIC

# The published-opioid sites the atlas recovered, keyed by the shortest observed peptide at
# each. These are the positive controls for M1. Mature f-numbers, matching docs/milk_atlas.md.
KNOWN_OPIOID_SITES = {
    "YPFPGPI":    "beta-casomorphin-7, bovine beta-casein f60",
    "YPFVEPIP":   "human beta-casomorphin-5 region, human beta-casein f51",
    "YPVEPFTES":  "neocasomorphin-6, bovine beta-casein f114",
    "YLGYLEQLL":  "alphas1-casein exorphin f91 (free-Tyr form)",
    "YGLFQISNK":  "alpha-lactorphin site, human alpha-lactalbumin f50",
}

# Aromatic-rich novel candidates: M2 is specifically about these.
AROMATIC_RICH = ("YPYAVWYYPQ", "YPTYHTHGRYVPP")

# Not a milk protein; excluded from every statistic and reported separately.
EXCLUDED_GENES = {"ACTA1", "ACTA2"}

PREREGISTERED = {
    "M1": "AUROC(known-opioid sites, novel sites) >= 0.70",
    "M2": "aromatic-rich novel candidates lose more score under exclusion volumes "
          "than the known controls do",
    "M3": "rho(score, aromatic count) weaker with exclusion than without "
          "(prior value without: +0.909)",
    "M4": "seed-to-seed Spearman rho >= 0.95 (S1 gate)",
}


MAX_EXACT_ORDERINGS = 200_000
PERMUTATION_SAMPLES = 20_000


def auroc_p(pos, neg, seed=11):
    """One-sided p for AUROC(pos, neg), exact when that is tractable, sampled otherwise.

    `D.exact_auroc_p` enumerates every assignment of which ranks are negatives, which is
    the honest statistic at dietary-set sizes: E8 has 3 antagonists, so C(17,3) = 680
    orderings. This set is balanced and much larger -- roughly 16 known against 17 novel,
    C(33,17) = 1.17e9 orderings, 1.7 million times more. Calling the exact function here
    does not return. (It was called here, and it did not return.)

    So: exact below MAX_EXACT_ORDERINGS, and above it a sampled permutation test, with
    which one was used recorded in the output rather than left to be inferred.
    """
    import random
    from math import comb
    npos, nneg = len(pos), len(neg)
    if not npos or not nneg:
        return {"p_one_sided": float("nan"), "method": "undefined (empty class)"}
    if comb(npos + nneg, nneg) <= MAX_EXACT_ORDERINGS:
        return {"p_one_sided": D.exact_auroc_p(pos, neg), "method": "exact"}
    observed = D.auroc(pos, neg)
    pool = list(pos) + list(neg)
    rng = random.Random(seed)
    hits = 0
    for _ in range(PERMUTATION_SAMPLES):
        rng.shuffle(pool)
        if D.auroc(pool[:npos], pool[npos:]) >= observed - 1e-9:
            hits += 1
    return {"p_one_sided": (hits + 1) / (PERMUTATION_SAMPLES + 1),
            "method": f"sampled permutation, n={PERMUTATION_SAMPLES}"}


def classify(seq):
    """known / novel, by whether the peptide starts with a published opioid peptide."""
    for ref in KNOWN_OPIOID_SITES:
        if seq.startswith(ref):
            return "known"
    return "novel"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default=Path("results/milk_atlas"), type=Path)
    ap.add_argument("--targets", default=Path("data/milk_atlas/milk_atlas_targets.tsv"),
                    type=Path)
    ap.add_argument("--out", default=Path("docs"), type=Path)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    with args.targets.open() as fh:
        meta = {r["sequence"]: r for r in csv.DictReader(fh, delimiter="\t")}

    runs = defaultdict(dict)
    for d in sorted(args.runs.glob("*_s*")):
        tag, _, seed = d.name.rpartition("_s")
        rows = D.load_run(d)
        if rows:
            # key by sequence, so excl and bare arms line up regardless of id
            runs[tag][seed] = {r["sequence"]: r for r in rows.values()
                               if r.get("sequence")}
    if not runs:
        raise SystemExit(f"no atlas runs under {args.runs}; "
                         f"run jobs/milk_atlas_screen.sbatch first")

    report = {
        "generated": date.today().isoformat(),
        "metric": METRIC,
        "design": "docs/milk_atlas.md",
        "preregistered": PREREGISTERED,
        "known_opioid_sites": KNOWN_OPIOID_SITES,
        "excluded_from_statistics": sorted(EXCLUDED_GENES),
        "configurations": {},
    }

    scores_by_tag = {}
    for tag, by_seed in sorted(runs.items()):
        seeds = sorted(by_seed)
        seqs = set.intersection(*(set(by_seed[s]) for s in seeds)) & set(meta)
        seqs = {s for s in seqs if meta[s]["gene"] not in EXCLUDED_GENES}
        if not seqs:
            continue
        score = {s: sum(float(by_seed[sd][s][METRIC]) for sd in seeds) / len(seeds)
                 for s in seqs}
        scores_by_tag[tag] = score

        # M4, seed stability
        rho_seeds = float("nan")
        if len(seeds) >= 2:
            shared = sorted(seqs)
            rho_seeds = D.spearman(
                [float(by_seed[seeds[0]][s][METRIC]) for s in shared],
                [float(by_seed[seeds[1]][s][METRIC]) for s in shared])

        known = [score[s] for s in seqs if classify(s) == "known"]
        novel = [score[s] for s in seqs if classify(s) == "novel"]
        arom = {s: D.n_aromatic(s[:8]) for s in seqs}

        m1 = None
        if known and novel:
            m1 = {
                "n_known": len(known), "n_novel": len(novel),
                "auroc": D.auroc(known, novel),
                **{f"p_{k}": v for k, v in auroc_p(known, novel).items()},
                "loo_known_auroc": sorted(
                    D.auroc([x for j, x in enumerate(known) if j != i], novel)
                    for i in range(len(known))) if len(known) > 1 else [],
                "median_known": D.median(known), "median_novel": D.median(novel),
            }

        report["configurations"][tag] = {
            "seeds": seeds,
            "n_scored": len(seqs),
            "M1_known_vs_novel": m1,
            "M3_rho_vs_aromatic": D.spearman([score[s] for s in sorted(seqs)],
                                             [arom[s] for s in sorted(seqs)]),
            "M4_seed_rho": rho_seeds,
            "M4_passes_S1_gate": bool(rho_seeds >= 0.95),
            "ranking": [
                {"rank": i, "sequence": s, "score": round(score[s], 4),
                 "class": classify(s), "known_as": KNOWN_OPIOID_SITES.get(
                     next((k for k in KNOWN_OPIOID_SITES if s.startswith(k)), ""), ""),
                 "gene": meta[s]["gene"], "species": meta[s].get("species", ""),
                 "site_mature": meta[s].get("site_position_mature", ""),
                 "n_aromatic_window8": arom[s]}
                for i, s in enumerate(sorted(seqs, key=lambda x: -score[x]), start=1)
            ],
        }

    # --- M2: exclusion-volume effect, per peptide -----------------------------------
    m2 = {}
    for excl_tag, bare_tag in (("a8excl", "a8bare"), ("b8excl", "b8bare")):
        if excl_tag not in scores_by_tag or bare_tag not in scores_by_tag:
            continue
        e, b = scores_by_tag[excl_tag], scores_by_tag[bare_tag]
        shared = sorted(set(e) & set(b))
        delta = {s: e[s] - b[s] for s in shared}
        known_d = [delta[s] for s in shared if classify(s) == "known"]
        arom_d = [delta[s] for s in shared if s in AROMATIC_RICH]
        m2[f"{excl_tag}_vs_{bare_tag}"] = {
            "median_delta_known": D.median(known_d) if known_d else None,
            "median_delta_novel": D.median([delta[s] for s in shared
                                            if classify(s) == "novel"]) or None,
            "delta_aromatic_rich": {s: round(delta[s], 4) for s in shared
                                    if s in AROMATIC_RICH},
            "aromatic_rich_lose_more_than_known": (
                bool(arom_d and known_d and D.median(arom_d) < D.median(known_d))),
            "per_peptide_delta": {s: round(delta[s], 4) for s in shared},
        }
    report["M2_exclusion_effect"] = m2

    out = args.out / "milk_atlas_result.json"
    out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")

    # --- console --------------------------------------------------------------------
    print(f"metric {METRIC}\n")
    for tag, cfg in report["configurations"].items():
        print(f"--- {tag}  seeds {','.join(cfg['seeds'])}  n={cfg['n_scored']} ---")
        m1 = cfg["M1_known_vs_novel"]
        if m1:
            loo = m1["loo_known_auroc"]
            loo_s = f"{min(loo):.3f}-{max(loo):.3f}" if loo else "n/a"
            verdict = "PASS" if m1["auroc"] >= 0.70 else "FAIL"
            print(f"  M1 known vs novel  n={m1['n_known']}/{m1['n_novel']}"
                  f"  AUROC {m1['auroc']:.3f} ({verdict} vs 0.70)"
                  f"  p={m1['p_p_one_sided']:.4f} ({m1['p_method']})  LOO {loo_s}")
        print(f"  M3 rho vs aromatic count   {cfg['M3_rho_vs_aromatic']:+.3f}")
        print(f"  M4 seed rho                {cfg['M4_seed_rho']:+.3f}"
              f"  {'PASS' if cfg['M4_passes_S1_gate'] else 'FAIL (<0.95) - no shortlist'}")
        print("  top 8:")
        for r in cfg["ranking"][:8]:
            label = r["known_as"] or "novel"
            print(f"    {r['rank']:>2}. {r['sequence'][:22]:24} {r['score']:.4f}"
                  f"  {r['class']:5} arom={r['n_aromatic_window8']}  {label[:44]}")
        print()
    for pair, v in report["M2_exclusion_effect"].items():
        print(f"M2 {pair}: median delta known={v['median_delta_known']}, "
              f"novel={v['median_delta_novel']}")
        print(f"   aromatic-rich deltas: {v['delta_aromatic_rich']}")
        print(f"   aromatic-rich lose more than known: "
              f"{v['aromatic_rich_lose_more_than_known']}")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
