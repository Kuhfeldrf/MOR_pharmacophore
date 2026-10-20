#!/usr/bin/env python3
"""Stage 1: build the pharmacophore template from the experimental complex.

Two sources, combined:

  PLIP     reads which receptor residues the peptide ligand actually engages, and by
           what kind of interaction, in the deposited pose.
  Align-it generates the ligand's own typed pharmacophore points at their 3D
           coordinates.

The template is Align-it's feature set; PLIP supplies the evidence for which features
are load-bearing, recorded alongside it. Decision D1 (PLAN.md): the template is read
from the experimental pose only, never a predicted one, because all three pose
generators tested in this project rotate the Tyr phenol into an H-bond the experimental
structure does not make.

Usage:
    python scripts/03_template.py --complex data/8f7q_complex.pdb --chain P \
        --out results/01_template
"""
import argparse
import json
import subprocess
import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from mor_pharmacophore.alignit import generate_pharmacophore, version  # noqa: E402
from mor_pharmacophore.plip_step import peptide_interactions  # noqa: E402


def protonate(in_path, out_sdf, ph=7.4):
    """Add hydrogens at `ph` and write SDF.

    Decision D4 (added 2026-10-20, after inspecting the first template). The deposited
    structure carries no hydrogens, so OpenBabel types the free N-terminal amine as a
    neutral donor and Align-it emits NO charge feature. The protonated ligand carries
    formal charge +1 and the template gains the POSC feature at the N-terminus.

    This matters more than a feature count: the N-terminal ammonium is the anchoring
    group the whole positional rule is built on, and a template without it would score
    candidates on everything except the interaction that defines the pharmacophore.
    """
    result = subprocess.run(
        ["obabel", str(in_path), "-O", str(out_sdf), "-p", str(ph)],
        capture_output=True, text=True)
    if result.returncode != 0 or not Path(out_sdf).exists():
        raise RuntimeError(f"obabel protonation failed: {result.stderr.strip()[:300]}")
    return Path(out_sdf)


def extract_chain(pdb_path, chain, out_path):
    """Write just `chain` of `pdb_path` as its own PDB, preserving coordinates."""
    kept = []
    for line in Path(pdb_path).read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            if len(line) > 21 and line[21] == chain:
                kept.append(line)
        elif line.startswith("TER") and kept:
            kept.append(line)
    if not kept:
        raise ValueError(f"no atoms on chain {chain} in {pdb_path}")
    Path(out_path).write_text("\n".join(kept) + "\nEND\n")
    return out_path


def parse_phar(phar_path):
    """Read an Align-it .phar file into feature dicts.

    Format is one feature per line: NAME X Y Z ALPHA [normal coords...], terminated by
    a line containing only '$$$$'.
    """
    features = []
    for line in Path(phar_path).read_text().splitlines():
        line = line.strip()
        if not line or line == "$$$$":
            continue
        parts = line.split()
        if len(parts) < 5:
            continue
        name = parts[0]
        try:
            x, y, z, alpha = (float(p) for p in parts[1:5])
        except ValueError:
            continue  # the first line is the molecule title
        features.append({"type": name, "x": x, "y": y, "z": z, "alpha": alpha})
    return features


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--complex", default=Path("data/8f7q_complex.pdb"), type=Path)
    ap.add_argument("--chain", default="P")
    ap.add_argument("--out", default=Path("results/01_template"), type=Path)
    ap.add_argument("--ph", default=7.4, type=float)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    print(f"align-it: {version()}")

    print(f"PLIP: reading interactions for chain {args.chain} ...")
    ix = peptide_interactions(args.complex, args.chain)

    ligand_pdb = extract_chain(args.complex, args.chain, args.out / "ligand.pdb")
    ligand_sdf = protonate(ligand_pdb, args.out / "ligand_h.sdf", ph=args.ph)
    print(f"ligand written to {ligand_pdb}, protonated at pH {args.ph} -> {ligand_sdf}")

    phar = args.out / "mor_template.phar"
    generate_pharmacophore(ligand_sdf, phar, mol_type="SDF")
    features = parse_phar(phar)
    print(f"template: {len(features)} features -> {phar}")

    counts = {}
    for f in features:
        counts[f["type"]] = counts.get(f["type"], 0) + 1

    record = {
        "generated": date.today().isoformat(),
        "alignit_version": version(),
        "source_complex": str(args.complex),
        "ligand_chain": args.chain,
        "protonation_pH": args.ph,
        "ligand": ix.ligand,
        "feature_counts": counts,
        "n_features": len(features),
        "features": features,
        "plip_interactions": {
            "n_interactions": ix.n_interactions,
            "hbond_residues": ix.hbond_residues,
            "hydrophobic_residues": ix.hydrophobic_residues,
            "saltbridge_residues": ix.saltbridge_residues,
            "pistack_residues": ix.pistack_residues,
        },
    }
    (args.out / "template.json").write_text(json.dumps(record, indent=2) + "\n")

    print(f"ligand: {ix.ligand}")
    print(f"feature counts: {counts}")
    if "POSC" not in counts:
        print("WARNING: no positive-charge feature in the template. The N-terminal "
              "ammonium anchor is missing; check protonation before trusting a score.")
    print(f"PLIP interactions: {ix.n_interactions}")
    print(f"  hbond       {ix.hbond_residues}")
    print(f"  hydrophobic {ix.hydrophobic_residues}")
    print(f"  saltbridge  {ix.saltbridge_residues}")
    print(f"  pistack     {ix.pistack_residues}")


if __name__ == "__main__":
    main()
