#!/usr/bin/env python3
"""Build the receptor-peptide complex from the deposited coordinates, with provenance.

Why this script exists
----------------------
`data/8f7q_complex.pdb` was inherited in the repository and described in PLAN.md as
"receptor chain R + chain P residues 1-5". That description was accurate about the FILE
and wrong about the STRUCTURE: the deposition 8F7Q resolves **21 residues** of chain P,

    YGGFM TSEKS QTPLV TLFKN A        (beta-endorphin 1-21)

Wang et al. (Cell 2023, 186:413-427), the paper that deposited it, describe exactly this:
residues M5-A21 are visible, with P13-A21 forming an alpha-helix turn at the TM1/2 cleft
and the hydrophilic M5-T12 stretch acting as a linker that "interact[s] mainly with
ECL2". So the template built from that file encoded the message motif alone and dropped
16 resolved ligand residues, including every contact with ECL2 - which the same paper
identifies as the element that "filter[s] specific opioid peptides".

The root cause was accepting an inherited coordinate file without checking it against the
deposition. Hence this script: the complex is now built from the deposited mmCIF, and its
provenance - source URL, residue inventory, and a SHA-256 of both input and output - is
written next to it so the same mistake cannot be made silently again.

Peptide extent is an explicit, recorded choice
----------------------------------------------
`--peptide-residues N` truncates the peptide to its first N residues, and N is written
into the provenance. Three extents are meaningful here:

  5   the message motif (YGGFM). What the old template was, by accident.
  8   message plus E8, whose side chain hydrogen-bonds R213(ECL2)/Q214(ECL2). This is
      the extent that MATCHES the 8-residue candidate window, so reference and candidate
      describe the same amount of molecule.
  21  everything resolved. Structurally complete, but a 21-residue reference cannot be
      covered by an 8-residue candidate window, and whole-peptide candidate conformers
      are not tractable (81 rotatable bonds for a 17-mer returns zero conformers).

The extent mismatch matters for the aromatic confound. With a 5-residue reference
carrying 2 AROM features, a candidate only has to land some aromatic ring on each of two
positions, and a candidate with more aromatic groups has more chances to do so - which is
the mechanism behind rho = +0.909 between score and aromatic count. Matching the
reference extent to the candidate window makes the comparison like-for-like.

Usage:
    python scripts/22_prepare_complex.py --cif data/8F7Q.cif \\
        --receptor-chain R --peptide-chain P --peptide-residues 8 \\
        --out data/8f7q_complex_r8.pdb
"""
import argparse
import hashlib
import json
from datetime import date
from pathlib import Path

# mmCIF atom_site fields used; auth_* are the numbering the PDB file and the paper use.
_WANT = ("group_PDB", "id", "type_symbol", "label_atom_id", "label_comp_id",
         "label_asym_id", "label_seq_id", "Cartn_x", "Cartn_y", "Cartn_z",
         "occupancy", "B_iso_or_equiv", "auth_seq_id", "auth_comp_id",
         "auth_asym_id", "auth_atom_id", "pdbx_PDB_model_num", "label_alt_id")


def parse_atom_site(cif_path):
    """Minimal mmCIF atom_site reader. Returns a list of dicts."""
    cols, rows, in_loop = [], [], False
    for line in Path(cif_path).read_text().splitlines():
        s = line.strip()
        if s.startswith("_atom_site."):
            cols.append(s.split(".", 1)[1].split()[0])
            in_loop = True
            continue
        if in_loop:
            if not s or s.startswith(("#", "loop_", "_")):
                break
            parts = s.split()
            if len(parts) < len(cols):
                continue
            rows.append(dict(zip(cols, parts)))
    return rows


def sha256(path):
    h = hashlib.sha256()
    h.update(Path(path).read_bytes())
    return h.hexdigest()


def pdb_line(i, a):
    """One ATOM/HETATM record in fixed-column PDB format."""
    name = a["atom"]
    # PDB atom-name column convention: 1-char elements start at column 14.
    nm = (" " + name).ljust(4) if len(name) < 4 and len(a["element"]) == 1 \
        else name.ljust(4)
    return ("%-6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s"
            % (a["group"], i, nm, " ", a["resname"][:3], a["chain"], a["resnum"],
               " ", a["x"], a["y"], a["z"], a["occ"], a["b"],
               a["element"].rjust(2)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cif", required=True, type=Path)
    ap.add_argument("--pdb-id", default=None)
    ap.add_argument("--source-url", default=None)
    ap.add_argument("--receptor-chain", default="R")
    ap.add_argument("--peptide-chain", default="P")
    ap.add_argument("--peptide-residues", default=0, type=int,
                    help="keep only the first N peptide residues; 0 = all resolved")
    ap.add_argument("--model", default="1")
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    rows = parse_atom_site(args.cif)
    if not rows:
        raise SystemExit(f"no atom_site records parsed from {args.cif}")

    def field(r, *names):
        for n in names:
            v = r.get(n)
            if v not in (None, ".", "?"):
                return v
        return None

    atoms = {args.receptor_chain: [], args.peptide_chain: []}
    for r in rows:
        if r.get("group_PDB") not in ("ATOM", "HETATM"):
            continue
        if field(r, "pdbx_PDB_model_num") not in (None, args.model):
            continue
        alt = field(r, "label_alt_id")
        if alt not in (None, "A"):           # keep one altloc only
            continue
        chain = field(r, "auth_asym_id", "label_asym_id")
        if chain not in atoms:
            continue
        try:
            resnum = int(field(r, "auth_seq_id", "label_seq_id"))
        except (TypeError, ValueError):
            continue
        atoms[chain].append({
            "group": r["group_PDB"],
            "atom": field(r, "auth_atom_id", "label_atom_id") or "",
            "resname": field(r, "auth_comp_id", "label_comp_id") or "UNK",
            "chain": chain, "resnum": resnum,
            "x": float(r["Cartn_x"]), "y": float(r["Cartn_y"]),
            "z": float(r["Cartn_z"]),
            "occ": float(field(r, "occupancy") or 1.0),
            "b": float(field(r, "B_iso_or_equiv") or 0.0),
            "element": (field(r, "type_symbol") or "C").strip(),
        })

    pep = atoms[args.peptide_chain]
    rec = atoms[args.receptor_chain]
    if not pep:
        raise SystemExit(f"no atoms for peptide chain {args.peptide_chain!r}")
    if not rec:
        raise SystemExit(f"no atoms for receptor chain {args.receptor_chain!r}")

    pep_res = sorted({(a["resnum"], a["resname"]) for a in pep})
    resolved = [(n, c) for n, c in pep_res]
    if args.peptide_residues:
        keep = {n for n, _ in resolved[:args.peptide_residues]}
        pep = [a for a in pep if a["resnum"] in keep]
        kept_res = [(n, c) for n, c in resolved if n in keep]
    else:
        kept_res = resolved

    args.out.parent.mkdir(parents=True, exist_ok=True)
    lines, i = [], 0
    for a in rec:
        i += 1
        lines.append(pdb_line(i, a))
    lines.append("TER")
    for a in pep:
        i += 1
        lines.append(pdb_line(i, a))
    lines.append("TER")
    lines.append("END")
    args.out.write_text("\n".join(lines) + "\n")

    three_to_one = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
        "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
        "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
        "TYR": "Y", "VAL": "V"}
    seq_resolved = "".join(three_to_one.get(c, "X") for _, c in resolved)
    seq_kept = "".join(three_to_one.get(c, "X") for _, c in kept_res)

    prov = {
        "generated": date.today().isoformat(),
        "pdb_id": args.pdb_id, "source_url": args.source_url,
        "source_cif": str(args.cif), "source_cif_sha256": sha256(args.cif),
        "model": args.model,
        "receptor_chain": args.receptor_chain,
        "n_receptor_residues": len({a["resnum"] for a in rec}),
        "n_receptor_atoms": len(rec),
        "peptide_chain": args.peptide_chain,
        "peptide_residues_resolved_in_deposition": len(resolved),
        "peptide_sequence_resolved": seq_resolved,
        "peptide_residues_kept": len(kept_res),
        "peptide_sequence_kept": seq_kept,
        "peptide_truncated": len(kept_res) < len(resolved),
        "output": str(args.out), "output_sha256": sha256(args.out),
        "note": ("Peptide extent is a recorded choice, not an accident. The previous "
                 "inherited file carried 5 of the 21 resolved peptide residues and that "
                 "was not checked against the deposition."),
    }
    pj = args.out.with_suffix(".provenance.json")
    pj.write_text(json.dumps(prov, indent=2) + "\n")

    print("wrote %s" % args.out)
    print("  receptor chain %s : %d residues, %d atoms"
          % (args.receptor_chain, prov["n_receptor_residues"], len(rec)))
    print("  peptide  chain %s : %d of %d resolved residues kept"
          % (args.peptide_chain, len(kept_res), len(resolved)))
    print("  resolved in deposition : %s" % seq_resolved)
    print("  kept in this file      : %s%s"
          % (seq_kept, "  (TRUNCATED)" if prov["peptide_truncated"] else ""))
    print("wrote %s" % pj)


if __name__ == "__main__":
    main()
