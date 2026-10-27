#!/usr/bin/env python3
"""What the dietary evaluation set can and cannot resolve.

Every validation statistic in this project is an agonist-vs-antagonist AUROC on the
dietary set, and the antagonist class has three members. This script measures what that
buys, on the project's own labels and scores, and it exists because the dual-template
result (`docs/dual_template_result.md`) turned on a comparison the set cannot make.

Four questions:

  Q1  Evidence floor. The exact permutation test enumerates every label assignment, so
      with 12 agonists and 3 antagonists there are C(15,3) = 455 of them and no result
      on this set can ever report a one-sided p below 1/455. The AUROC is likewise
      granular at 1/(12*3), so one discordant pair moves it by a full 0.028.

  Q2  Can a working template be detected at all? Power of the one-sided permutation test
      to reject AUROC = 0.5, as a function of the antagonist count and the true AUROC.

  Q3  Is the observed A-vs-B gap distinguishable from zero? Paired bootstrap over
      peptides, resampling agonists and antagonists independently, on the real scores.

  Q4  Can ANY two templates be adjudicated here? Power of a paired bootstrap to resolve
      a true AUROC gap between two templates whose score vectors correlate at rho - the
      regime this project is actually in, since T3 measured rho = 0.917.

Q4 is the one that bears on design decisions: if it is near zero at n_ant = 3, then no
template swap can be validated on this set and the evaluation set is the binding
constraint, not the template.

Usage:
    python scripts/28_power_analysis.py --runs results/dietary \\
        --template-a p8excl --template-b b8excl --out docs
"""
import argparse
import csv
import itertools
import json
from datetime import date
from math import comb, sqrt
from pathlib import Path

import numpy as np
from scipy.stats import norm

SEED = 20261027


def auroc(pos, neg):
    """P(pos > neg) with ties at 0.5 - the Mann-Whitney statistic."""
    pos, neg = np.asarray(pos)[:, None], np.asarray(neg)[None, :]
    return float((pos > neg).mean() + 0.5 * (pos == neg).mean())


def exact_p(pos, neg):
    """Exact one-sided permutation p, enumerating every label assignment."""
    alls = np.concatenate([pos, neg])
    n, k = len(alls), len(neg)
    obs = auroc(pos, neg)
    hits = tot = 0
    for combo in itertools.combinations(range(n), k):
        mask = np.zeros(n, bool)
        mask[list(combo)] = True
        hits += auroc(alls[~mask], alls[mask]) >= obs - 1e-12
        tot += 1
    return hits / tot, tot


def hanley_se(auc, n_pos, n_neg):
    """Hanley & McNeil (1982) standard error of an AUC.

    Variance goes as 1/n_pos + 1/n_neg, so the SMALLER class dominates it. With three
    antagonists, adding agonists cannot buy precision - which is the single most
    important consequence of this set's composition.
    """
    q1 = auc / (2 - auc)
    q2 = 2 * auc ** 2 / (1 + auc)
    return sqrt((auc * (1 - auc)
                 + (n_pos - 1) * (q1 - auc ** 2)
                 + (n_neg - 1) * (q2 - auc ** 2)) / (n_pos * n_neg))


def logit_ci(auc, se, z=1.96):
    """95% CI computed on the logit scale and back-transformed.

    At this n the symmetric Wald interval runs past 1.0 - the failure Nicholls (2014,
    JCAMD 28:887, 10.1007/s10822-014-9753-z) works through explicitly, prescribing this
    transform as the fix. The project's percentile bootstrap (23_dietary_screen.py)
    cannot exceed 1.0 but is correspondingly optimistic when resamples separate
    perfectly, so both are reported.
    """
    if not 0 < auc < 1:
        return [float("nan"), float("nan")]
    l = np.log(auc / (1 - auc))
    sl = se / (auc * (1 - auc))
    return [float(1 / (1 + np.exp(-(l + s * z * sl)))) for s in (-1, 1)]


def read_seed_averaged(runs, tag):
    """Mean tversky_ref_mean per peptide across every <tag>_s<seed> run present."""
    per = {}
    for d in sorted(Path(runs).glob(f"{tag}_s*")):
        tsv = d / "block0000.tsv"
        if not tsv.exists():
            continue
        with tsv.open() as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                per.setdefault(r["id"], []).append(float(r["tversky_ref_mean"]))
    if not per:
        raise SystemExit(f"no runs for tag {tag!r} under {runs}")
    return {k: float(np.mean(v)) for k, v in per.items()}


def power_single(rng, n_ago, n_ant, true_auc, alpha=0.05, reps=600, perms=600):
    """Power to reject AUROC = 0.5 one-sided at alpha."""
    delta = sqrt(2) * norm.ppf(true_auc)
    wins = 0
    for _ in range(reps):
        g = rng.normal(delta, 1.0, n_ago)
        t = rng.normal(0.0, 1.0, n_ant)
        obs = auroc(g, t)
        pool = np.concatenate([g, t])
        null = np.array([auroc(p[:n_ago], p[n_ago:])
                         for p in (rng.permutation(pool) for _ in range(perms))])
        wins += (1 + (null >= obs - 1e-12).sum()) / (perms + 1) <= alpha
    return wins / reps


def power_paired(rng, n_ago, n_ant, auc_a, auc_b, rho, reps=400, boots=400):
    """Power of a paired bootstrap to resolve a true gap between correlated templates."""
    da, db = sqrt(2) * norm.ppf(auc_a), sqrt(2) * norm.ppf(auc_b)
    wins = 0
    for _ in range(reps):
        z1, z2 = rng.normal(size=n_ago), rng.normal(size=n_ago)
        ga = da + z1
        gb = db + rho * z1 + sqrt(1 - rho ** 2) * z2
        w1, w2 = rng.normal(size=n_ant), rng.normal(size=n_ant)
        ta = w1
        tb = rho * w1 + sqrt(1 - rho ** 2) * w2
        bs = np.empty(boots)
        for j in range(boots):
            gi = rng.integers(0, n_ago, n_ago)
            ti = rng.integers(0, n_ant, n_ant)
            bs[j] = auroc(ga[gi], ta[ti]) - auroc(gb[gi], tb[ti])
        lo, hi = np.percentile(bs, [2.5, 97.5])
        wins += (lo > 0) or (hi < 0)
    return wins / reps


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default=Path("results/dietary"), type=Path)
    ap.add_argument("--meta", default=Path("data/dietary/dietary_peptides.tsv"), type=Path)
    ap.add_argument("--template-a", default="p8excl")
    ap.add_argument("--template-b", default="b8excl")
    ap.add_argument("--out", default=Path("docs"), type=Path)
    args = ap.parse_args()

    rng = np.random.default_rng(SEED)
    with args.meta.open() as fh:
        meta = {r["id"]: r for r in csv.DictReader(fh, delimiter="\t")}

    # The primary pool of docs/dual_template_design.md 6: mu-dominant assays, published,
    # supported termini. Reconstructed here rather than imported so this script can be
    # read on its own.
    pool = [i for i, m in meta.items()
            if m["assay_receptor_bias"] == "mu_dominant"
            and m["published"] == "yes" and m["terminal_supported"] == "yes"]
    ago = [i for i in pool if meta[i]["activity"] == "Agonist"]
    ant = [i for i in pool if meta[i]["activity"] == "Antagonist"]

    sa = read_seed_averaged(args.runs, args.template_a)
    sb = read_seed_averaged(args.runs, args.template_b)
    auc_a = auroc([sa[i] for i in ago], [sa[i] for i in ant])
    auc_b = auroc([sb[i] for i in ago], [sb[i] for i in ant])
    pa, tot = exact_p(np.array([sa[i] for i in ago]), np.array([sa[i] for i in ant]))
    pb, _ = exact_p(np.array([sb[i] for i in ago]), np.array([sb[i] for i in ant]))

    report = {
        "generated": str(date.today()),
        "seed": SEED,
        "pool": {"n_agonist": len(ago), "n_antagonist": len(ant)},
        "Q1_evidence_floor": {
            f"{na}x{nt}": {"label_assignments": comb(na + nt, nt),
                           "min_one_sided_p": 1 / comb(na + nt, nt),
                           "auroc_granule": 1 / (na * nt)}
            for na, nt in ((12, 3), (14, 3), (12, 5), (12, 8), (12, 12))
        },
        "observed": {
            "auroc_a": auc_a, "exact_p_a": pa,
            "auroc_b": auc_b, "exact_p_b": pb,
            "delta": auc_a - auc_b,
            "discordant_pairs": round((auc_a - auc_b) * len(ago) * len(ant)),
            "label_assignments": tot,
            "hanley_se_a": hanley_se(auc_a, len(ago), len(ant)),
            "hanley_se_b": hanley_se(auc_b, len(ago), len(ant)),
            "logit_ci95_a": logit_ci(auc_a, hanley_se(auc_a, len(ago), len(ant))),
            "logit_ci95_b": logit_ci(auc_b, hanley_se(auc_b, len(ago), len(ant))),
        },
    }
    for nm, au in (("a", auc_a), ("b", auc_b)):
        se = hanley_se(au, len(ago), len(ant))
        lo, hi = logit_ci(au, se)
        print(f"   {nm.upper()}8 AUROC {au:.4f}  Hanley SE {se:.3f}  "
              f"logit CI95 [{lo:.3f}, {hi:.3f}]")

    print(f"pool: {len(ago)} agonists / {len(ant)} antagonists")
    print(f"Q1 evidence floor: min one-sided p = {1/tot:.5f} over {tot} assignments; "
          f"AUROC granule {1/(len(ago)*len(ant)):.4f}")
    print(f"   observed A {auc_a:.4f} (p {pa:.5f})   B {auc_b:.4f} (p {pb:.5f})   "
          f"delta {auc_a-auc_b:+.4f}")

    print("Q2 power to reject AUROC=0.5 (12 agonists):")
    report["Q2_power_single"] = {}
    for n_ant in (3, 5, 8, 12):
        row = {f"{a}": power_single(rng, 12, n_ant, a) for a in (0.75, 0.85, 0.95)}
        report["Q2_power_single"][f"n_ant={n_ant}"] = row
        print(f"   n_ant {n_ant:2d}  " + "  ".join(f"true {k}: {v:.2f}" for k, v in row.items()))

    # Q3 on the real scores: resample agonists and antagonists independently.
    aa = np.array([sa[i] for i in ago]); ta = np.array([sa[i] for i in ant])
    ab = np.array([sb[i] for i in ago]); tb = np.array([sb[i] for i in ant])
    d = np.empty(20000)
    for i in range(d.size):
        gi = rng.integers(0, len(ago), len(ago))
        ti = rng.integers(0, len(ant), len(ant))
        d[i] = auroc(aa[gi], ta[ti]) - auroc(ab[gi], tb[ti])
    lo, hi = (float(x) for x in np.percentile(d, [2.5, 97.5]))
    report["Q3_delta_bootstrap"] = {
        "delta": auc_a - auc_b, "ci95": [lo, hi],
        "p_delta_le_0": float((d <= 0).mean()), "reps": int(d.size),
    }
    print(f"Q3 delta {auc_a-auc_b:+.4f}  CI95 [{lo:+.4f}, {hi:+.4f}]  "
          f"P(delta<=0) {(d <= 0).mean():.3f}")

    print("Q4 power to adjudicate two templates correlated at rho=0.9 (12 agonists):")
    report["Q4_power_paired_rho0.9"] = {}
    for n_ant in (3, 5, 8, 12, 20):
        row = {f"gap={g}": power_paired(rng, 12, n_ant, 0.90, 0.90 - g, 0.9)
               for g in (0.05, 0.10, 0.15)}
        report["Q4_power_paired_rho0.9"][f"n_ant={n_ant}"] = row
        print(f"   n_ant {n_ant:2d}  " + "  ".join(f"{k}: {v:.2f}" for k, v in row.items()))

    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "power_analysis.json").write_text(json.dumps(report, indent=2))
    print(f"wrote {args.out / 'power_analysis.json'}")


if __name__ == "__main__":
    main()
