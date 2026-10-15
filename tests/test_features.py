import pytest

from mor_pharmacophore.features import feature_counts, feature_positions

MET_ENK = "YGGFM"       # met-enkephalin
BCM7 = "YPFPGPI"        # beta-casomorphin-7
NEGATIVE = "AAAAA"      # no aromatic side chains


def test_opioid_peptides_share_two_aromatic_features():
    assert feature_counts(MET_ENK)["Aromatic"] == 2   # Tyr1, Phe4
    assert feature_counts(BCM7)["Aromatic"] == 2      # Tyr1, Phe3


def test_negative_control_has_no_aromatic_feature():
    assert feature_counts(NEGATIVE)["Aromatic"] == 0


def test_positive_charge_from_free_nterminus():
    assert feature_counts(MET_ENK)["PosIonizable"] >= 1


def test_feature_positions_are_3d_points():
    pts = feature_positions(MET_ENK, "Aromatic")
    assert len(pts) == 2 and all(len(p) == 3 for p in pts)
