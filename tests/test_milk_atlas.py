import importlib.util
from pathlib import Path

import pytest

_SPEC = importlib.util.spec_from_file_location(
    "milk_atlas", Path(__file__).resolve().parents[1] / "scripts" / "25_milk_motif_atlas.py")
atlas = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(atlas)


def make_index(proteins):
    ix = atlas.ProteomeIndex()
    for acc, seq in proteins:
        ix.add(acc, f"{acc}_TEST", acc, "test protein", seq, "test", "Testus testus")
    return ix.build()


def test_locate_reports_protein_and_1based_offset():
    ix = make_index([("A", "MKWYPFPGPIQQ"), ("B", "AAAAYGGFLAAA")])
    (prot, start), = ix.locate("YPFPGPI")
    assert prot["acc"] == "A" and start == 3          # 0-based; 4th residue
    (prot, start), = ix.locate("YGGFL")
    assert prot["acc"] == "B" and start == 4


def test_peptide_spanning_the_separator_is_not_a_hit():
    """The join separator must never let a match straddle two proteins."""
    ix = make_index([("A", "MKQQY"), ("B", "PFPGPIQQ")])
    assert ix.locate("YPFPGPI") == []
    # and the separator itself is not matchable as sequence
    assert ix.locate(atlas.ProteomeIndex.SEP) == []


def test_locate_finds_every_occurrence_in_one_protein():
    ix = make_index([("A", "YPFPAAAYPFPAAA")])
    hits = ix.locate("YPFP")
    assert [s for _, s in hits] == [0, 7]


def test_tier1_motifs_match_the_known_opioid_peptides():
    """Each tier-1 class must match the reference peptides that motivate it."""
    import re
    byclass = {label: pattern for label, pattern, _ in atlas.TIER1}
    assert re.match(byclass["message_GGF"], "YGGFM")
    assert re.match(byclass["message_GGF"], "FGGFTGARK")      # N/OFQ, Phe anchor
    assert re.match(byclass["casomorphin_YP"], "YPFPGPI")     # beta-casomorphin-7
    assert re.match(byclass["casomorphin_YP"], "YPWF")        # endomorphin-1
    assert re.match(byclass["casomorphin_YP"], "YPVEPF")      # neocasomorphin-6
    assert re.match(byclass["casoxin_YPxY"], "YPSYGLN")       # casoxin A
    assert re.match(byclass["lactorphin_YxxF"], "YGLF")       # alpha-lactorphin
    assert re.match(byclass["lactorphin_YxxF"], "YLLF")       # beta-lactorphin


def test_tier1_does_not_match_a_non_opioid_start():
    import re
    for label, pattern, _ in atlas.TIER1:
        assert not re.match(pattern, "AAAAAAA"), label
        assert not re.match(pattern, "GGFYM"), label          # no N-terminal aromatic anchor


def test_scan_sites_reports_overlapping_matches():
    """Lookahead scanning must not consume the match, or adjacent sites are lost."""
    sites = atlas.scan_sites(make_index([("A", "YPYPYPF")]), atlas.TIER1, "1_canonical")
    positions = sorted(s["site_position_1based"] for s in sites
                       if s["motif_class"] == "casomorphin_YP")
    assert len(positions) >= 2 and positions[0] == 1


def test_mature_position_applies_the_signal_peptide():
    # beta-casomorphin-7 sits at precursor 75 of P02666 and is published as f60-66.
    assert atlas.mature_position("P02666", 75) == (60, "validated")
    # alpha-lactorphin site: precursor 69 of P00709 -> f50.
    assert atlas.mature_position("P00709", 69)[0] == 50


def test_mature_position_strips_the_isoform_suffix():
    assert atlas.mature_position("P47710-3", 137)[0] == atlas.mature_position("P47710", 137)[0]


def test_unknown_accession_gets_no_mature_position_rather_than_a_guess():
    assert atlas.mature_position("P68133", 71) == ("", "")


def test_site_inside_the_signal_peptide_has_no_mature_coordinate():
    assert atlas.mature_position("P02666", 5)[0] == ""


def test_species_is_parsed_and_distinguishes_the_lactalbumins():
    _, _, _, _, human = atlas.parse_header(
        "sp|P00709|LALBA_HUMAN Alpha-lactalbumin OS=Homo sapiens GN=LALBA PE=1 SV=1")
    _, _, _, _, bovine = atlas.parse_header(
        "sp|P00711|LALBA_BOVIN Alpha-lactalbumin OS=Bos taurus GN=LALBA PE=1 SV=2")
    assert human == "Homo sapiens" and bovine == "Bos taurus"


def test_parse_header_pulls_accession_and_gene():
    acc, name, gene, desc, _ = atlas.parse_header(
        "sp|P02666|CASB_BOVIN Beta-casein OS=Bos taurus GN=CSN2 PE=1 SV=2")
    assert (acc, name, gene) == ("P02666", "CASB_BOVIN", "CSN2")
    assert desc == "Beta-casein"
