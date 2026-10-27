#!/usr/bin/env python3
"""Which ligand residues touch which receptor regions, straight from the coordinates.

Written to settle a documented claim rather than to explore: `docs/structural_basis.md` 5
asserted that the selectivity-determining region of beta-endorphin sits at residues
~13-21 and is therefore out of reach of any tractable template. This script measures the
contacts and reports which ligand residue reaches ECL2, the TM6 tip, ECL1 and the TM1
tip, so that claim can be checked instead of inherited. See `docs/post_t3_design.md` 2.

Deliberately stdlib-only, so it can be run with `python3 -I` against a file downloaded
from the PDB without that file's directory being on the import path.

Regions are given in the numbering of the structure being read. 8F7Q and 8F7R use human
P35372 numbering (D149 = 3.32); 6DDE/6DDF/9WSW and most other MOR entries use the
398-residue convention (D147 = 3.32), a constant -2 offset. Passing the wrong region map
is the obvious way to get a plausible, wrong answer, so --regions is explicit.

Usage:
    python3 -I scripts/29_ligand_contacts.py 8f7q.cif --ligand P --receptor R \\
        --out docs/8f7q_contacts.tsv
"""
import argparse
import gzip
from collections import defaultdict
from math import dist

# MOR regions, 8F7Q/human P35372 numbering.
REGIONS_8F7Q = [
    ("TM1-tip", range(60, 80)),
    ("ECL1", range(126, 140)),
    ("pocket-core", list(range(143, 156)) + list(range(236, 245)) + list(range(290, 298))),
    ("ECL2", range(205, 225)),
    ("TM5-tip", range(228, 236)),
    ("TM6-tip", range(298, 310)),
    ("TM7/pocket", range(316, 335)),
]
REGION_MAPS = {
    "8f7q": REGIONS_8F7Q,
    # The same regions shifted by -2 for the 398-residue convention.
    "mor398": [(n, range(r.start - 2, r.stop - 2)) if isinstance(r, range)
               else (n, [x - 2 for x in r]) for n, r in REGIONS_8F7Q],
}


def region_of(num, regions):
    for name, rng in regions:
        if num in rng:
            return name
    return "other"


def load_atoms(path):
    """Heavy atoms of model 1 from an mmCIF, as dicts. Column order is read, not assumed."""
    op = gzip.open if str(path).endswith(".gz") else open
    cols, rows, hdr = [], [], False
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith("loop_"):
                cols, hdr = [], False
                continue
            if line.startswith("_atom_site."):
                cols.append(line.strip().split(".", 1)[1])
                hdr = True
                continue
            if hdr and (line.startswith("ATOM") or line.startswith("HETATM")):
                rows.append(line.split())
            elif hdr and rows:
                break
    ix = {c: i for i, c in enumerate(cols)}
    need = ("auth_asym_id", "auth_seq_id", "label_comp_id", "label_atom_id",
            "Cartn_x", "Cartn_y", "Cartn_z", "type_symbol", "pdbx_PDB_model_num")
    missing = [c for c in need if c not in ix]
    if missing:
        raise SystemExit(f"{path}: _atom_site is missing {missing}")
    out = []
    for r in rows:
        if len(r) < len(cols) or r[ix["pdbx_PDB_model_num"]] != "1":
            continue
        if r[ix["type_symbol"]] == "H":
            continue
        try:
            xyz = (float(r[ix["Cartn_x"]]), float(r[ix["Cartn_y"]]), float(r[ix["Cartn_z"]]))
        except ValueError:
            continue
        out.append({"chain": r[ix["auth_asym_id"]], "resi": r[ix["auth_seq_id"]],
                    "resn": r[ix["label_comp_id"]], "xyz": xyz})
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("structure", help="mmCIF path (.cif or .cif.gz)")
    ap.add_argument("--ligand", required=True, help="auth chain of the peptide ligand")
    ap.add_argument("--receptor", required=True, help="auth chain of the receptor")
    ap.add_argument("--cutoff", type=float, default=4.5)
    ap.add_argument("--regions", default="8f7q", choices=sorted(REGION_MAPS))
    ap.add_argument("--out", default=None, help="TSV to write; stdout only if omitted")
    args = ap.parse_args()

    regions = REGION_MAPS[args.regions]
    atoms = load_atoms(args.structure)
    lig = [a for a in atoms if a["chain"] == args.ligand]
    rec = [a for a in atoms if a["chain"] == args.receptor]
    if not lig or not rec:
        raise SystemExit(f"empty selection: ligand {len(lig)} atoms, receptor {len(rec)}")

    closest = defaultdict(dict)
    for la in lig:
        for ra in rec:
            d = dist(la["xyz"], ra["xyz"])
            if d <= args.cutoff:
                key = (int(ra["resi"]), ra["resn"])
                if key not in closest[int(la["resi"])] or d < closest[int(la["resi"])][key]:
                    closest[int(la["resi"])][key] = d

    resn_of = {int(a["resi"]): a["resn"] for a in lig}
    reach = defaultdict(set)
    lines = ["lig_resi\tlig_resn\tregion\trec_resn\trec_resi\tmin_dist_A"]
    for resi in sorted(resn_of):
        hits = closest.get(resi, {})
        if not hits:
            print(f"  lig {resi:3d} {resn_of[resi]:4s}  -- no contact --")
            lines.append(f"{resi}\t{resn_of[resi]}\t-\t-\t-\t-")
            continue
        parts = []
        for (rnum, rnm), d in sorted(hits.items(), key=lambda kv: kv[1]):
            reg = region_of(rnum, regions)
            reach[reg].add(resi)
            lines.append(f"{resi}\t{resn_of[resi]}\t{reg}\t{rnm}\t{rnum}\t{d:.2f}")
            parts.append(f"[{reg}] {rnm}{rnum} {d:.1f}")
        print(f"  lig {resi:3d} {resn_of[resi]:4s}  " + "; ".join(parts[:5]))

    print("\n=== ligand residues reaching each region ===")
    for name, _ in regions:
        if name in reach:
            print(f"  {name:12s} <- {sorted(reach[name])}")
    if "other" in reach:
        print(f"  {'other':12s} <- {sorted(reach['other'])}")

    if args.out:
        with open(args.out, "w") as fh:
            fh.write("\n".join(lines) + "\n")
        print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
