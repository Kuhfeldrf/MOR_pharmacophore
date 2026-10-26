#!/usr/bin/env python3
"""Opioid-motif atlas over the milk proteome, intersected with observed peptides.

Two questions, and they are not the same question:

  1. **Where are the opioid motifs in the milk proteins?** A sequence-level scan of the
     parent proteome. Every site found here is a *potential* opioid peptide N-terminus.
  2. **Which of those are actually released, with the motif at the peptide's N-terminus?**
     This needs peptidomics, not sequence. `data/library/peptides.tsv` is 14,571 peptides
     observed by MS in infant gastric and intestinal digesta, so it answers it directly.

Question 2 is the one that matters, because of the positional rule
(`mor_pharmacophore.nterm_rule`): opioid binding needs a FREE alpha-amine at the anchoring
Tyr. A peptide that merely *contains* `YPFP` has that amine locked in a peptide bond and
cannot anchor. So a motif site is only interesting if digestion actually cuts immediately
before it. This script reports the three populations separately:

    released_nterminal   observed peptide starts exactly at the motif site  -> candidate
    spanning             observed peptide contains the site, starts upstream -> blocked
    not_observed         site never seen in any observed peptide             -> no evidence

Motif definitions are derived from this repository's own reference peptides, not invented
here. Tier 1 is the set of N-terminal patterns actually exhibited by known opioid peptides
in `data/dietary/dietary_peptides.tsv` and `data/reference/wang2023_tableS2.csv`. Tier 2 is
the mechanistic generalisation the source paper supports: `docs/structural_basis.md` 4
quotes Wang et al. finding that YGGF (endorphin/dynorphin), YPW (endomorphin), Y[D-Ala]F
(deltorphin) and FGGF (N/OFQ) all occupy one "highly aromatic pocket", i.e. what the
receptor reads is an aromatic anchor plus a second aromatic two or three residues later.
Tier 2 is therefore a superset and is reported, but Tier 1 is what gets screened.

Caveat carried in the output, not just here: Tier 2 is permissive by design and a Tier-2
hit is a hypothesis, not a prediction. Nothing in this script measures binding; it selects
what is worth scoring.

Usage:
    python scripts/25_milk_motif_atlas.py \\
        --fasta human=/path/HumanMilkProteinDatabase_v2.fasta \\
        --fasta bovine=/path/CowMilkProteinDatabase_v2.fasta \\
        --library data/library/peptides.tsv --out data/milk_atlas
"""
import argparse
import json
import re
from collections import defaultdict
from datetime import date
from pathlib import Path

# Tier 1: N-terminal patterns exhibited by known opioid peptides in this repo.
# Each entry is (label, regex, the reference peptides that motivate it).
TIER1 = [
    ("message_GGF", r"[YF]GGF",
     "Met-/Leu-enkephalin YGGFM/YGGFL, beta-endorphin YGGFMTSEK..., "
     "dynorphin A YGGFLRRIRPKLK, N/OFQ FGGFTGARK... (F anchor)"),
    ("casomorphin_YP", r"YP[FWYVLIM]",
     "bovine beta-casomorphins YPFP/YPFPG/YPFPGPI, human YPFV, soymorphin YPFVVN, "
     "neocasomorphin YPVEPF, endomorphin-1 YPWF, casoxin B YPYY"),
    ("casoxin_YPxY", r"YP[STAG][YFW]",
     "casoxin A YPSYGLN, casoxin-6 SRYPSY (internal site)"),
    ("lactorphin_YxxF", r"Y[GL][LGSF][FY]",
     "alpha-lactorphin YGLF, beta-lactorphin YLLF, lactoferroxin A YLGSGY"),
]

# Tier 2: the aromatic-chamber generalisation. Aromatic anchor, second aromatic at +2 or
# +3. Deliberately permissive; a superset of Tier 1.
TIER2 = [
    ("aromatic_pair", r"[YF].{1,2}[FWY]",
     "docs/structural_basis.md 4: YGGF, YPW, Y[D-Ala]F and FGGF share one aromatic pocket"),
]

# The milk proteins an opioid peptide could plausibly come from in quantity. Sites outside
# these are still reported (the scan is proteome-wide) but flagged, because a motif in a
# trace signalling protein is not a dietary opioid source in any practical sense.
MAJOR_GENES = {
    "CSN1S1", "CSN1S2", "CSN2", "CSN3",        # caseins
    "LALBA", "PAEP", "BLG", "LGB",             # alpha-lactalbumin, beta-lactoglobulin
    "LTF", "ALB", "MUC1", "OPN", "SPP1",       # lactoferrin, serum albumin, mucin, osteopontin
    "CSN1S2A", "CSN1S2B", "XDH", "BTN", "FOLR1",
}

# Signal-peptide length per accession, for converting precursor numbering to the MATURE
# numbering the literature uses.
#
# This is not cosmetic. UniProt sequences are precursors including the signal peptide, so
# a site this script reports at position 75 of P02666 is the peptide every paper calls
# f60-66 -- beta-casomorphin-7. Reporting only precursor coordinates would put every
# position in this atlas off by the signal-peptide length relative to the literature,
# silently, in the direction that makes a known peptide look novel.
#
# `validated` entries are confirmed by this repo: the known opioid peptide at that site
# lands on its published mature position once the offset is applied (see docs/milk_atlas.md
# 4). `uniprot` entries are taken from the UniProt SIGNAL annotation and are not
# independently checked here. An accession absent from this table gets a blank mature
# position rather than a guess.
SIGNAL_PEPTIDE = {
    # accession: (length, basis)
    "P02666": (15, "validated"),   # CASB_BOVIN    YPFPGPI -> f60-66, beta-casomorphin-7
    "P02662": (15, "validated"),   # CASA1_BOVIN   YLGYLE  -> f91-96, alphas1 exorphin
    "P00709": (19, "validated"),   # LALBA_HUMAN   YGLF    -> f50-53, alpha-lactorphin site
    "P00711": (19, "uniprot"),     # LALBA_BOVIN
    "P05814": (15, "uniprot"),     # CASB_HUMAN
    "P47710": (15, "uniprot"),     # CASA1_HUMAN
    "P02663": (15, "uniprot"),     # CASA2_BOVIN
    "P02668": (21, "uniprot"),     # CASK_BOVIN
    "P07498": (20, "uniprot"),     # CASK_HUMAN
    "P02754": (16, "uniprot"),     # LACB_BOVIN, beta-lactoglobulin
    "P24627": (19, "uniprot"),     # TRFL_BOVIN, lactoferrin
    "P02788": (19, "uniprot"),     # TRFL_HUMAN, lactoferrin
}


def mature_position(accession, precursor_pos):
    """(mature_position, basis) for a precursor position, or ("", "") if unknown.

    Isoform suffixes are stripped: P47710-3 uses P47710's signal peptide.
    """
    base = accession.split("-")[0]
    entry = SIGNAL_PEPTIDE.get(base)
    if not entry:
        return "", ""
    length, basis = entry
    pos = precursor_pos - length
    # A site inside the signal peptide has no mature coordinate, and a signal-peptide
    # fragment is not a dietary peptide anyway.
    return (pos, basis) if pos >= 1 else ("", basis)


def read_fasta(path):
    """Yield (header, sequence). Minimal reader: these are plain UniProt FASTA files."""
    header, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:], []
            elif line:
                chunks.append(line)
    if header is not None:
        yield header, "".join(chunks)


def parse_header(h):
    """Pull accession, entry name and gene out of a UniProt-style header."""
    acc, name = "", h.split()[0] if h else ""
    parts = h.split("|")
    if len(parts) >= 3:
        acc, name = parts[1], parts[2].split()[0]
    gene = ""
    m = re.search(r"GN=([A-Za-z0-9_.-]+)", h)
    if m:
        gene = m.group(1)
    desc = h.split(" ", 1)[1] if " " in h else ""
    desc = re.split(r" (?:OS|OX|GN|PE|SV)=", desc)[0]
    # Species from OS=, falling back to the entry-name suffix (LALBA_HUMAN). Worth
    # carrying explicitly: P00709 is LALBA_HUMAN and P00711 is LALBA_BOVIN, so an atlas
    # that reports only the gene name cannot distinguish a breast-milk peptide from a
    # formula one -- which is the whole point of the comparison this feeds.
    species = ""
    m = re.search(r"OS=([^=]+?)(?: [A-Z]{2}=|$)", h)
    if m:
        species = m.group(1).strip()
    elif "_" in name:
        species = name.rsplit("_", 1)[1]
    return acc, name, gene, desc, species


class ProteomeIndex:
    """Concatenated proteome with offset mapping, for fast peptide location.

    14,571 peptides against ~900 proteins is 13 million naive substring searches. Joining
    the proteome into one string and using `str.find` pushes the search into C and makes
    the whole mapping a few seconds. The separator must be a character that cannot occur
    in a sequence, or a peptide could be "found" straddling two proteins.
    """

    SEP = "|"

    def __init__(self):
        self.proteins = []          # (acc, name, gene, desc, seq, offset, source)
        self._parts, self._cursor = [], 0

    def add(self, acc, name, gene, desc, seq, source, species=""):
        self.proteins.append({
            "acc": acc, "name": name, "gene": gene, "desc": desc,
            "species": species, "seq": seq, "offset": self._cursor,
            "source": source})
        self._parts.append(seq)
        self._cursor += len(seq) + len(self.SEP)

    def build(self):
        self.joined = self.SEP.join(self._parts)
        self._starts = [p["offset"] for p in self.proteins]
        return self

    def locate(self, peptide):
        """All (protein, start_index_0based) where `peptide` occurs."""
        hits, pos = [], self.joined.find(peptide)
        while pos != -1:
            # binary search for the protein containing `pos`
            lo, hi = 0, len(self._starts) - 1
            while lo < hi:
                mid = (lo + hi + 1) // 2
                if self._starts[mid] <= pos:
                    lo = mid
                else:
                    hi = mid - 1
            prot = self.proteins[lo]
            start = pos - prot["offset"]
            if 0 <= start and start + len(peptide) <= len(prot["seq"]):
                hits.append((prot, start))
            pos = self.joined.find(peptide, pos + 1)
        return hits


def scan_sites(index, motifs, tier):
    """Every motif occurrence in every protein."""
    sites = []
    for prot in index.proteins:
        for label, pattern, basis in motifs:
            for m in re.finditer(f"(?=({pattern}))", prot["seq"]):
                start = m.start()
                matched = m.group(1)
                mat, mat_basis = mature_position(prot["acc"], start + 1)
                sites.append({
                    "source": prot["source"], "accession": prot["acc"],
                    "entry": prot["name"], "gene": prot["gene"],
                    "species": prot["species"], "protein": prot["desc"],
                    "major_milk_protein": "yes" if prot["gene"] in MAJOR_GENES else "no",
                    "tier": tier, "motif_class": label, "motif_regex": pattern,
                    "site_position_1based": start + 1,
                    "site_position_mature": mat,
                    "mature_numbering_basis": mat_basis,
                    "matched": matched,
                    "context_8mer": prot["seq"][start:start + 8],
                    "basis": basis,
                })
    return sites


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", action="append", required=True, metavar="LABEL=PATH",
                    help="repeatable; e.g. human=/path/HumanMilk.fasta")
    ap.add_argument("--library", default=Path("data/library/peptides.tsv"), type=Path)
    ap.add_argument("--out", default=Path("data/milk_atlas"), type=Path)
    ap.add_argument("--min-length", default=4, type=int,
                    help="shortest peptide to carry forward as a screening candidate; "
                         "the template needs an anchor plus a second aromatic, so a "
                         "tripeptide cannot express the motif")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    index = ProteomeIndex()
    for spec in args.fasta:
        if "=" not in spec:
            raise SystemExit(f"--fasta expects LABEL=PATH, got {spec!r}")
        label, path = spec.split("=", 1)
        n = 0
        for header, seq in read_fasta(path):
            if not seq:
                continue
            acc, name, gene, desc, species = parse_header(header)
            index.add(acc, name, gene, desc, seq, label, species)
            n += 1
        print(f"{label}: {n} proteins from {path}")
    index.build()

    sites = scan_sites(index, TIER1, "1_canonical") + scan_sites(index, TIER2, "2_aromatic_pair")
    print(f"motif sites: {len(sites)} "
          f"(tier 1: {sum(1 for s in sites if s['tier'].startswith('1'))})")

    # ---- observed peptides -> parent protein and position --------------------------
    import csv
    with args.library.open() as fh:
        library = list(csv.DictReader(fh, delimiter="\t"))
    print(f"library: {len(library)} observed peptides")

    # site lookup: (accession, position) -> list of site records
    by_site = defaultdict(list)
    for s in sites:
        by_site[(s["accession"], s["site_position_1based"])].append(s)

    released, spanning = [], []
    unmapped = 0
    for row in library:
        pep = row["sequence"]
        hits = index.locate(pep)
        if not hits:
            unmapped += 1
            continue
        for prot, start in hits:
            start1 = start + 1
            # Released at a motif site: the peptide's own N-terminus IS the site, so the
            # anchoring Tyr carries a free alpha-amine (nterm_rule.supports_site).
            for s in by_site.get((prot["acc"], start1), []):
                released.append({**s, "peptide_id": row["peptide_id"],
                                 "peptide": pep, "peptide_length": len(pep),
                                 "peptide_start_1based": start1,
                                 "matrices": row.get("matrices", ""),
                                 "phases": row.get("phases", ""),
                                 "relation": "released_nterminal"})
            # Spanning: contains a site but starts upstream of it -> amine blocked.
            for off in range(1, len(pep)):
                for s in by_site.get((prot["acc"], start1 + off), []):
                    if off + len(s["matched"]) <= len(pep):
                        spanning.append({**s, "peptide_id": row["peptide_id"],
                                         "peptide": pep,
                                         "peptide_start_1based": start1,
                                         "offset_into_peptide": off,
                                         "relation": "spanning_blocked_nterm"})

    def write_tsv(path, rows, fields=None):
        if not rows:
            path.write_text("")
            return 0
        fields = fields or list(rows[0].keys())
        with path.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                               extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)
        return len(rows)

    write_tsv(args.out / "motif_sites.tsv", sites)
    write_tsv(args.out / "released_nterminal.tsv", released)
    write_tsv(args.out / "spanning_blocked.tsv", spanning)

    # ---- screening target table: Tier 1, released, long enough, deduplicated -------
    seen, targets = set(), []
    for r in sorted(released, key=lambda r: (r["tier"], r["gene"], r["site_position_1based"])):
        if not r["tier"].startswith("1"):
            continue
        if r["peptide_length"] < args.min_length:
            continue
        if r["peptide"] in seen:
            continue
        seen.add(r["peptide"])
        targets.append({
            "id": f"milk__{r['gene'] or r['accession']}__{r['site_position_1based']}__{r['peptide_id']}",
            "sequence": r["peptide"],
            "set": "milk_atlas",
            "c_term": "free",
            "motif_class": r["motif_class"],
            "gene": r["gene"],
            "accession": r["accession"],
            "site_position_1based": r["site_position_1based"],
            "site_position_mature": r["site_position_mature"],
            "species": r["species"],
            "major_milk_protein": r["major_milk_protein"],
            "source": r["source"],
        })
    write_tsv(args.out / "milk_atlas_targets.tsv", targets,
              fields=["id", "sequence", "set", "c_term", "motif_class", "gene",
                      "accession", "species", "site_position_1based",
                      "site_position_mature", "major_milk_protein", "source"])

    # ---- summary -------------------------------------------------------------------
    def tally(rows, key):
        out = defaultdict(int)
        for r in rows:
            out[r[key]] += 1
        return dict(sorted(out.items(), key=lambda kv: -kv[1]))

    t1_sites = [s for s in sites if s["tier"].startswith("1")]
    t1_released = [r for r in released if r["tier"].startswith("1")]
    summary = {
        "generated": date.today().isoformat(),
        "n_proteins": len(index.proteins),
        "n_observed_peptides": len(library),
        "n_observed_unmapped_to_proteome": unmapped,
        "motif_definitions": {
            "tier1_canonical": [{"class": a, "regex": b, "basis": c} for a, b, c in TIER1],
            "tier2_aromatic_pair": [{"class": a, "regex": b, "basis": c} for a, b, c in TIER2],
        },
        "sites": {
            "total": len(sites),
            "tier1": len(t1_sites),
            "tier1_in_major_milk_proteins": sum(
                1 for s in t1_sites if s["major_milk_protein"] == "yes"),
            "tier1_by_class": tally(t1_sites, "motif_class"),
            "tier1_by_gene_top20": dict(list(tally(t1_sites, "gene").items())[:20]),
        },
        "released_nterminal": {
            "rows": len(released),
            "tier1_rows": len(t1_released),
            "tier1_distinct_peptides": len({r["peptide"] for r in t1_released}),
            "tier1_distinct_sites": len({(r["accession"], r["site_position_1based"])
                                         for r in t1_released}),
            "tier1_by_class": tally(t1_released, "motif_class"),
            "tier1_by_gene": tally(t1_released, "gene"),
        },
        "spanning_blocked": {
            "rows": len(spanning),
            "distinct_peptides": len({r["peptide"] for r in spanning}),
        },
        "screening_targets": {
            "n": len(targets),
            "rule": "tier 1 motif, released as peptide N-terminus, "
                    f"length >= {args.min_length}, deduplicated by sequence",
            "by_class": tally(targets, "motif_class"),
            "by_gene": tally(targets, "gene"),
        },
        "caveats": [
            "Tier 2 is permissive by construction (aromatic anchor + second aromatic at "
            "+2/+3) and a tier-2 hit is a hypothesis, not a prediction. Only tier 1 "
            "becomes a screening target.",
            "A motif site is evidence of a potential N-terminus, not of release. The "
            "released_nterminal set is the only one with peptidomic evidence behind it.",
            "Nothing here measures binding. The atlas selects what is worth scoring.",
            "Peptides are mapped to parents by exact substring match, so a peptide shared "
            "between two parents is reported against both.",
        ],
    }
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    print(f"\ntier 1 sites: {len(t1_sites)} "
          f"({summary['sites']['tier1_in_major_milk_proteins']} in major milk proteins)")
    print(f"tier 1 by class: {summary['sites']['tier1_by_class']}")
    print(f"released at a tier-1 site: {summary['released_nterminal']['tier1_distinct_peptides']} "
          f"distinct peptides at {summary['released_nterminal']['tier1_distinct_sites']} sites")
    print(f"  by gene: {summary['released_nterminal']['tier1_by_gene']}")
    print(f"spanning but blocked: {summary['spanning_blocked']['distinct_peptides']} distinct peptides")
    print(f"screening targets written: {len(targets)} -> {args.out/'milk_atlas_targets.tsv'}")
    print(f"unmapped observed peptides: {unmapped}")


if __name__ == "__main__":
    main()
