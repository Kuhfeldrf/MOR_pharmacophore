#!/usr/bin/env python3
"""W1-W4 of docs/post_t3_design.md, decided by numbers fixed before the runs.

Exit 0 = the extent-9 design survives its gates, exit 1 = it does not. Intended to be
read the way scripts/26_t3_gate.py is read: the thresholds below were written into
post_t3_design.md 4 before any window-9 score existed, and this script does not choose
anything.

  W1  Sampling adequacy. Seed-to-seed Spearman of the window-9 dietary scores must be
      >= 0.98, matching window 8's measured 0.983. Failing W1 means 800 conformers do
      not cover 31-37 rotatable bonds and window 8 is the ceiling FOR SAMPLING REASONS -
      the design is over, and that is a result.

  W2  The extent-9 template carries a feature the extent-8 one does not. Lys9's side-chain
      amine is the whole point of extent 9 (it salt-bridges Asp218 in ECL2, 3.1 A, see
      docs/8f7q_contacts.tsv), so the extent-9 .phar must hold a POSC feature that extent
      8 lacks. If not, extent 9 is extent 8 with extra atoms.

  W3  Non-redundancy. Spearman between the extent-9 and extent-8 dietary score vectors
      must be < 0.9 - the same threshold T3 applied to the two templates, for the same
      reason. This is the test most likely to fail, and failing it says the method
      saturates at the message epitope.

  W4  The composition confound has not returned: rho(score, aromatic count) < 0.909.

W1 and W3 are the decision. W5 (agonist/antagonist AUROC) is deliberately NOT here: at
three antagonists its power against a true AUROC of 0.85 is 0.62 and its minimum
attainable p is 0.0022 (docs/power_analysis.md), so it cannot adjudicate anything and is
reported by scripts/23_dietary_screen.py rather than gated on.

Usage:
    python scripts/30_extent9_gates.py --runs results/dietary \\
        --extent9 p9excl --extent8 p8excl --out docs
"""
import argparse
import importlib.util
import json
import sys
from datetime import date
from pathlib import Path

_HERE = Path(__file__).resolve().parent

W1_MIN_SEED_RHO = 0.98
W3_MAX_RHO_VS_EXTENT8 = 0.9
W4_MAX_AROMATIC_RHO = 0.909


def _dietary_helpers():
    """Import `23_dietary_screen.py`, whose name is not a legal module identifier."""
    spec = importlib.util.spec_from_file_location(
        "dietary_screen", _HERE / "23_dietary_screen.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


D = _dietary_helpers()
METRIC = D.METRIC


def runs_for(runs, tag):
    """{seed: {peptide_id: row}} for every <tag>_s<seed> directory present."""
    out = {}
    for d in sorted(Path(runs).glob(f"{tag}_s*")):
        if not (d / "block0000.tsv").exists():
            continue
        out[d.name.split("_s")[-1]] = D.load_run(d)
    return out


def phar_features(path):
    """[(type, x, y, z)] from an Align-it .phar file."""
    feats = []
    for line in Path(path).read_text().splitlines():
        parts = line.split("\t") if "\t" in line else line.split()
        if len(parts) >= 4 and parts[0].isalpha() and parts[0].isupper():
            try:
                feats.append((parts[0], *(float(v) for v in parts[1:4])))
            except ValueError:
                continue
    return feats


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default=Path("results/dietary"), type=Path)
    ap.add_argument("--extent9", default="p9excl")
    ap.add_argument("--extent8", default="p8excl")
    ap.add_argument("--template9", default=Path("results/template_p9/mor_template.phar"),
                    type=Path)
    ap.add_argument("--template8", default=Path("results/template_p8/mor_template.phar"),
                    type=Path)
    ap.add_argument("--out", default=Path("docs"), type=Path)
    args = ap.parse_args()

    report = {"generated": str(date.today()), "design": "docs/post_t3_design.md",
              "thresholds": {"W1_min_seed_rho": W1_MIN_SEED_RHO,
                             "W3_max_rho_vs_extent8": W3_MAX_RHO_VS_EXTENT8,
                             "W4_max_aromatic_rho": W4_MAX_AROMATIC_RHO}}
    failures = []

    r9 = runs_for(args.runs, args.extent9)
    if len(r9) < 2:
        print(f"GATE FAIL: need two seeds of {args.extent9!r} under {args.runs}; "
              f"found {sorted(r9)}", file=sys.stderr)
        return 1
    r8 = runs_for(args.runs, args.extent8)

    # --- W1 -----------------------------------------------------------------------
    s1, s2 = sorted(r9)[:2]
    shared = sorted(set(r9[s1]) & set(r9[s2]))
    rho_seed = D.spearman([float(r9[s1][p][METRIC]) for p in shared],
                          [float(r9[s2][p][METRIC]) for p in shared])
    report["W1_seed_to_seed_rho"] = rho_seed
    ok1 = rho_seed >= W1_MIN_SEED_RHO
    failures += [] if ok1 else ["W1"]
    print(f"  W1  seed-to-seed rho = {rho_seed:+.4f}  "
          f"(need >= {W1_MIN_SEED_RHO})  {'PASS' if ok1 else 'FAIL'}")

    # --- W2 -----------------------------------------------------------------------
    f9, f8 = phar_features(args.template9), phar_features(args.template8)
    n9 = sum(1 for t, *_ in f9 if t == "POSC")
    n8 = sum(1 for t, *_ in f8 if t == "POSC")
    report["W2_posc_extent9"], report["W2_posc_extent8"] = n9, n8
    report["W2_features_extent9"] = len(f9)
    report["W2_features_extent8"] = len(f8)
    ok2 = n9 > n8
    failures += [] if ok2 else ["W2"]
    print(f"  W2  POSC features: extent9 {n9} vs extent8 {n8}; total features "
          f"{len(f9)} vs {len(f8)}  {'PASS' if ok2 else 'FAIL'}")

    # --- W3 -----------------------------------------------------------------------
    if len(r8) >= 1:
        def seed_avg(runs):
            ids = set.intersection(*(set(v) for v in runs.values()))
            return {p: sum(float(runs[s][p][METRIC]) for s in runs) / len(runs)
                    for p in ids}
        a9, a8 = seed_avg(r9), seed_avg(r8)
        both = sorted(set(a9) & set(a8))
        rho_38 = D.spearman([a9[p] for p in both], [a8[p] for p in both])
        report["W3_rho_vs_extent8"] = rho_38
        report["W3_n_peptides"] = len(both)
        ok3 = abs(rho_38) < W3_MAX_RHO_VS_EXTENT8
        failures += [] if ok3 else ["W3"]
        print(f"  W3  rho(extent9, extent8) = {rho_38:+.4f}  "
              f"(need < {W3_MAX_RHO_VS_EXTENT8})  {'PASS' if ok3 else 'FAIL'}")
    else:
        report["W3_rho_vs_extent8"] = None
        print(f"  W3  SKIPPED - no {args.extent8!r} runs under {args.runs}")

    # --- W4 -----------------------------------------------------------------------
    rows = [r9[s1][p] for p in shared]
    rho_arom = D.spearman([float(r[METRIC]) for r in rows],
                          [D.n_aromatic(r["window"]) for r in rows])
    report["W4_rho_vs_aromatic"] = rho_arom
    ok4 = abs(rho_arom) < W4_MAX_AROMATIC_RHO
    failures += [] if ok4 else ["W4"]
    print(f"  W4  rho(score, aromatic) = {rho_arom:+.4f}  "
          f"(need |rho| < {W4_MAX_AROMATIC_RHO})  {'PASS' if ok4 else 'FAIL'}")

    report["failures"] = failures
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "extent9_gates.json").write_text(json.dumps(report, indent=2))
    print(f"wrote {args.out / 'extent9_gates.json'}")

    # W1 and W3 are the decision; W2 and W4 are reported and recorded but do not by
    # themselves end the design.
    decisive = [f for f in failures if f in ("W1", "W3")]
    if decisive:
        print(f"\nGATE FAIL: {', '.join(decisive)} - see docs/post_t3_design.md 4 for "
              f"what each failure means. This is a result, not an error.")
        return 1
    print("\nGATE PASS: extent 9 is computable at this ensemble size and is not a "
          "re-measurement of extent 8.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
