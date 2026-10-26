#!/usr/bin/env python3
"""Decide whether the full-library screen under template B is worth running.

Exit 0 = run it, exit 1 = do not. Intended as the condition in a job chain, so the
decision is made by a number rather than by whoever is at the terminal.

The gate is T3 from `docs/dual_template_design.md` 6. The full library is 14,595 targets
at 800 conformers, roughly 52 CPU-hours per seed, and its only purpose is the E4'
percentile (T5). If the two templates' dietary score vectors correlate at rho >= 0.9 then
template B is re-measuring what template A already measured: the dual-template design is
unjustified, the consensus is not a consensus of anything, and that compute buys a
percentile against a distribution nobody should be quoting. The cheap screen comes first
precisely so this can be checked for the price of a few minutes.

A missing or unreadable result file exits 1 (do not run), not 0. The expensive default is
never the fallback.

Usage:
    python scripts/26_t3_gate.py --result docs/dual_template_result.json [--threshold 0.9]
"""
import argparse
import json
import sys
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--result", default=Path("docs/dual_template_result.json"), type=Path)
    ap.add_argument("--threshold", default=0.9, type=float,
                    help="rho at or above which the templates count as redundant")
    args = ap.parse_args()

    if not args.result.exists():
        print(f"GATE FAIL: {args.result} does not exist - run the dietary screen and "
              f"scripts/24_dual_template_consensus.py first", file=sys.stderr)
        return 1
    try:
        report = json.loads(args.result.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        print(f"GATE FAIL: cannot read {args.result}: {exc}", file=sys.stderr)
        return 1

    blocks = {"seed_averaged": report.get("seed_averaged")}
    blocks.update({f"seed {s}": b for s, b in (report.get("per_seed") or {}).items()})
    rhos = {k: v.get("T3_rho_between_templates")
            for k, v in blocks.items() if isinstance(v, dict)}
    rhos = {k: v for k, v in rhos.items() if isinstance(v, (int, float))}
    if not rhos:
        print(f"GATE FAIL: no T3_rho_between_templates in {args.result}", file=sys.stderr)
        return 1

    for k, v in rhos.items():
        print(f"  T3 rho({k}) = {v:+.3f}")

    # Judge on the seed-averaged value when present; it is the primary statistic.
    primary = rhos.get("seed_averaged", max(rhos.values(), key=abs))
    if abs(primary) >= args.threshold:
        print(f"GATE FAIL: |rho| = {abs(primary):.3f} >= {args.threshold} - the two "
              f"templates are redundant, so the ~52 CPU-hour library screen would buy "
              f"nothing. Revisit the design (docs/dual_template_design.md 9) instead.")
        return 1
    print(f"GATE PASS: |rho| = {abs(primary):.3f} < {args.threshold} - templates carry "
          f"different information; the library screen is justified.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
