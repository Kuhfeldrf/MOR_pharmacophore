#!/usr/bin/env python3
"""Add receptor-derived exclusion volumes to the pharmacophore template (D2.4).

Why
---
The template built by `03_template.py` has 16 features, all attractive, and no exclusion
volumes. Align-it reserves score-table column 6 for pharmacophore/exclusion-sphere
overlap, and in every run of this pipeline so far that column is identically zero:
nothing penalises a candidate for bulk, or for features the pocket cannot physically
accommodate. That is the mechanistic origin of the aromatic-count confound - a His-rich
8-mer accumulates template coverage with no counter-pressure, because the score only ever
asks what the candidate *covers*, never what it would *clash with*.

Two independent literature lines say this is the missing piece:

  - Both relevant precedents include receptor-derived exclusion volumes. The 4DKL MOR
    tetrapeptide screen used LigandScout's defaults (PMC6865014); apo2ph4 keeps all of
    them, removing only duplicates (PMC9832483).

  - Horvath, Marcou & Varnek (J Chem Inf Model 2013, 53:1543, "Do Not Hesitate to Use
    Tversky...") distinguish `growable` queries, which admit active analogues carrying
    extra features and are best served by query-biased Tversky (alpha > 0.5), from
    `conservative` queries that allow no growth, whose specific application domain is
    *candidate*-biased Tversky (alpha < 0.5). The opioid message pocket is sterically
    tight, so it is a conservative query - and this pipeline has been scoring it with a
    strongly query-biased metric. Exclusion volumes are the 3D-native way to encode
    `no growth allowed`, rather than inverting the metric asymmetry.

Note the Horvath caveat, which matters for how much weight that second argument carries:
their descriptors are ISIDA colored fragment / pharmacophore *triplet counts*, i.e.
topological 2D, not 3D overlap. The growable-vs-conservative distinction is a conceptual
import here, not a transferred benchmark result.

How
---
Exclusion spheres are placed on receptor heavy atoms that line the ligand pocket:
within `--contact` A of any ligand heavy atom, excluding the ligand's own chain. Spheres
closer together than `--dedupe` A are merged to one representative, following apo2ph4's
practice of keeping all exclusion volumes but removing duplicates.

Align-it's `alpha` is a Gaussian width, not a hard radius, so it is swept rather than
guessed: a larger alpha is a tighter sphere. Templates are written for each requested
alpha so the effect can be measured instead of assumed.

Scoring must then pass `--withExclusion`, or the spheres are parsed and ignored.

Prediction, recorded before any scoring (see DECISIONS.md D2.4)
--------------------------------------------------------------
Exclusion volumes should penalise over-large, feature-dense windows more than
message-like ones: E6 (score vs aromatic count) should fall while E1 (the potency pair)
holds. If both fall together, the collinearity found in `16_reselect_with_e6.py` is
intrinsic to the method rather than an artefact of an incomplete template. Either result
is informative; the point of writing it down now is that neither can be claimed as the
expected one afterwards.

Usage:
    python scripts/20_template_exclusion.py --complex data/8f7q_complex.pdb \\
        --ligand-chain P --template results/01_template/mor_template.phar \\
        --out results/01_template --alphas 0.5 1.0
"""
import argparse
import json
from datetime import date
from pathlib import Path

LIGAND_PHAR_CODES = ("AROM", "HDON", "HACC", "LIPO", "CHARGE", "POSC", "NEGC",
                     "HYBH", "HYBL")


def read_pdb_atoms(path):
    """(chain, resname, resnum, atomname, element, x, y, z) for heavy atoms."""
    out = []
    for line in Path(path).read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        name = line[12:16].strip()
        element = (line[76:78].strip() or name[:1]).upper()
        if element == "H":
            continue
        try:
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
        except ValueError:
            continue
        out.append((line[21], line[17:20].strip(), line[22:26].strip(), name,
                    element, x, y, z))
    return out


def read_phar(path):
    """Existing template feature lines, unchanged, plus a parsed code list."""
    lines, codes = [], []
    for raw in Path(path).read_text().splitlines():
        s = raw.strip()
        if not s or s == "$$$$":
            continue
        if "\t" not in s and " " not in s:
            lines.append(raw)                      # the leading name line
            continue
        code = s.split()[0]
        if code in LIGAND_PHAR_CODES:
            codes.append(code)
        lines.append(raw)
    return lines, codes


def d2(a, b):
    return (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--complex", default=Path("data/8f7q_complex.pdb"), type=Path)
    ap.add_argument("--ligand-chain", default="P")
    ap.add_argument("--template",
                    default=Path("results/01_template/mor_template.phar"), type=Path)
    ap.add_argument("--out", default=Path("results/01_template"), type=Path)
    ap.add_argument("--contact", default=5.0, type=float,
                    help="receptor heavy atom within this distance of any ligand "
                         "heavy atom is pocket-lining")
    ap.add_argument("--dedupe", default=1.5, type=float,
                    help="merge exclusion spheres closer than this")
    ap.add_argument("--clear", default=2.0, type=float,
                    help="drop any sphere this close to an existing ATTRACTIVE "
                         "template feature, so exclusion volumes cannot sit on top of "
                         "the interactions the template is built from")
    ap.add_argument("--alphas", default=[0.5, 1.0], type=float, nargs="+")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    atoms = read_pdb_atoms(args.complex)
    lig = [a for a in atoms if a[0] == args.ligand_chain]
    rec = [a for a in atoms if a[0] != args.ligand_chain]
    if not lig:
        raise SystemExit(f"no ligand atoms in chain {args.ligand_chain!r}")
    if not rec:
        raise SystemExit("no receptor atoms found outside the ligand chain")

    lig_xyz = [(a[5], a[6], a[7]) for a in lig]
    base_lines, base_codes = read_phar(args.template)

    # attractive feature coordinates, to keep exclusion spheres off them
    feat_xyz = []
    for raw in base_lines:
        p = raw.split()
        if p and p[0] in LIGAND_PHAR_CODES and len(p) >= 4:
            feat_xyz.append((float(p[1]), float(p[2]), float(p[3])))

    contact2, clear2, dedupe2 = args.contact ** 2, args.clear ** 2, args.dedupe ** 2

    # ---- pocket-lining receptor atoms --------------------------------------------
    lining = []
    for a in rec:
        p = (a[5], a[6], a[7])
        if any(d2(p, q) <= contact2 for q in lig_xyz):
            lining.append(a)

    # ---- drop those sitting on an attractive feature ------------------------------
    kept = [a for a in lining
            if not any(d2((a[5], a[6], a[7]), f) <= clear2 for f in feat_xyz)]
    dropped_on_feature = len(lining) - len(kept)

    # ---- dedupe --------------------------------------------------------------------
    spheres = []
    for a in kept:
        p = (a[5], a[6], a[7])
        if any(d2(p, s[:3]) < dedupe2 for s in spheres):
            continue
        spheres.append((p[0], p[1], p[2], a[1], a[2], a[3]))

    residues = sorted({f"{s[3]}{s[4]}" for s in spheres})

    written = []
    for alpha in args.alphas:
        tag = ("%g" % alpha).replace(".", "p")
        path = args.out / f"mor_template_excl_a{tag}.phar"
        lines = list(base_lines)
        for s in spheres:
            # Align-it EXCL: code x y z alpha hasNormal nx ny nz (no normal -> zeros)
            lines.append("EXCL\t%.3f\t%.3f\t%.3f\t%g\t0\t0\t0\t0"
                         % (s[0], s[1], s[2], alpha))
        path.write_text("\n".join(lines) + "\n$$$$\n")
        written.append(str(path))
        print("wrote %s  (%d attractive + %d EXCL, alpha=%g)"
              % (path, len(base_codes), len(spheres), alpha))

    meta = {
        "generated": date.today().isoformat(),
        "source_complex": str(args.complex), "ligand_chain": args.ligand_chain,
        "base_template": str(args.template),
        "base_feature_counts": {c: base_codes.count(c) for c in sorted(set(base_codes))},
        "parameters": {"contact_A": args.contact, "dedupe_A": args.dedupe,
                       "clear_of_feature_A": args.clear, "alphas": args.alphas},
        "n_receptor_heavy_atoms": len(rec),
        "n_pocket_lining": len(lining),
        "n_dropped_on_attractive_feature": dropped_on_feature,
        "n_exclusion_spheres": len(spheres),
        "lining_residues": residues,
        "templates_written": written,
        "scoring_requirement": "align_it must be called with --withExclusion",
    }
    (args.out / "template_exclusion.json").write_text(json.dumps(meta, indent=2) + "\n")

    print("\nreceptor heavy atoms      %d" % len(rec))
    print("pocket-lining (<= %.1f A)   %d" % (args.contact, len(lining)))
    print("dropped on a feature      %d" % dropped_on_feature)
    print("exclusion spheres kept    %d" % len(spheres))
    print("lining residues (%d): %s" % (len(residues), ", ".join(residues)))
    print("\nNOTE scoring must pass --withExclusion or the spheres are ignored.")


if __name__ == "__main__":
    main()
