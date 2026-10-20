#!/usr/bin/env python3
"""Stage 3: ingest the digestion peptide library into the screening library of record.

Reads the Dallas lab combined peptidomics workbook and writes one row per unique
peptide, carrying the sheet(s) it came from and the digestive phase it was observed in.
Phase is carried from here rather than rejoined later, because the narrative prioritises
peptides present in the same digestive phase as the receptor.

Nothing is filtered on opioid likeness. The library is the whole library; selection
happens downstream where it can be measured.

Usage:
    python scripts/01_library.py --xlsx <workbook.xlsx> --out data/library
"""
import argparse
import hashlib
import json
import re
from collections import defaultdict
from datetime import date
from pathlib import Path

import openpyxl

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
MIN_LENGTH = 2

# Sheet name -> digestive phase. The abbott sheet is per-column instead (see below).
SHEET_PHASE = {
    "Preterm_gastric": "gastric",
    "pre-term_intestinal": "intestinal",
    "term_intestinial": "intestinal",
}

# abbott intensity columns -> (matrix, phase). A non-empty cell means the peptide was
# observed in that condition.
ABBOTT_COLUMNS = {
    "Human Milk Intestinal": ("human_milk", "intestinal"),
    "Human Milk Gastric": ("human_milk", "gastric"),
    "Human Milk Undigested": ("human_milk", "undigested"),
    "Bovine Undigested": ("bovine", "undigested"),
    "Bovine Gastric": ("bovine", "gastric"),
    "Bovine Intestinal": ("bovine", "intestinal"),
}


def peptide_id(sequence):
    """Stable id from the sequence itself, so ids do not depend on row order."""
    return "pep_" + hashlib.sha256(sequence.encode()).hexdigest()[:10]


def is_standard(sequence):
    return bool(sequence) and set(sequence) <= STANDARD_AA


def read_workbook(xlsx_path):
    """Return (records, rejects, sheet_counts).

    records: sequence -> {sheets, matrices, phases}
    rejects: list of (sheet, raw value, reason) - recorded, never silently dropped.
    """
    wb = openpyxl.load_workbook(xlsx_path, read_only=True, data_only=True)
    records = defaultdict(lambda: {"sheets": set(), "matrices": set(), "phases": set()})
    rejects = []
    sheet_counts = {}

    for ws in wb.worksheets:
        rows = ws.iter_rows(values_only=True)
        header = next(rows, None)
        if header is None:
            sheet_counts[ws.title] = 0
            continue
        header = [str(h).strip() if h is not None else "" for h in header]
        kept = 0

        for row in rows:
            if not row or row[0] is None:
                continue
            raw = str(row[0]).strip().upper()
            if not is_standard(raw):
                rejects.append((ws.title, raw[:50], "non-standard residue"))
                continue
            if len(raw) < MIN_LENGTH:
                rejects.append((ws.title, raw, f"shorter than {MIN_LENGTH} aa"))
                continue

            rec = records[raw]
            rec["sheets"].add(ws.title)
            kept += 1

            if ws.title == "abbott":
                for col_name, (matrix, phase) in ABBOTT_COLUMNS.items():
                    if col_name not in header:
                        continue
                    value = row[header.index(col_name)]
                    if value is not None and value != "":
                        rec["matrices"].add(matrix)
                        rec["phases"].add(phase)
            else:
                phase = SHEET_PHASE.get(ws.title)
                if phase:
                    rec["phases"].add(phase)
                # The preterm/term sheets are human milk digestion datasets.
                rec["matrices"].add("human_milk")

        sheet_counts[ws.title] = kept

    return records, rejects, sheet_counts


def write_library(records, out_dir):
    out_dir.mkdir(parents=True, exist_ok=True)
    tsv = out_dir / "peptides.tsv"
    columns = ["peptide_id", "sequence", "length", "starts_with_Y",
               "source_sheets", "matrices", "phases"]

    with tsv.open("w") as fh:
        fh.write("\t".join(columns) + "\n")
        for sequence in sorted(records):
            rec = records[sequence]
            fh.write("\t".join([
                peptide_id(sequence),
                sequence,
                str(len(sequence)),
                "yes" if sequence.startswith("Y") else "no",
                ";".join(sorted(rec["sheets"])),
                ";".join(sorted(rec["matrices"])) or "NA",
                ";".join(sorted(rec["phases"])) or "NA",
            ]) + "\n")
    return tsv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--xlsx", required=True, type=Path)
    ap.add_argument("--out", default=Path("data/library"), type=Path)
    args = ap.parse_args()

    records, rejects, sheet_counts = read_workbook(args.xlsx)
    tsv = write_library(records, args.out)

    tyr = sum(1 for s in records if s.startswith("Y"))
    provenance = {
        "generated": date.today().isoformat(),
        "source_workbook": args.xlsx.name,
        "source_sha256": hashlib.sha256(args.xlsx.read_bytes()).hexdigest(),
        "sheet_counts": sheet_counts,
        "unique_peptides": len(records),
        "tyr_starting": tyr,
        "rejected": len(rejects),
        "rejects": [{"sheet": s, "value": v, "reason": r} for s, v, r in rejects],
        "min_length": MIN_LENGTH,
    }
    (args.out / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")

    print(f"wrote {tsv}")
    for sheet, n in sheet_counts.items():
        print(f"  {sheet:24s} {n:7d} rows kept")
    print(f"  {'unique peptides':24s} {len(records):7d}")
    print(f"  {'Tyr-starting':24s} {tyr:7d}")
    print(f"  {'rejected':24s} {len(rejects):7d}")


if __name__ == "__main__":
    main()
