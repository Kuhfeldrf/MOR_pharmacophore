"""The provenance block written into each screen's block JSON.

A seed is not provenance if the library stack that consumed it is unrecorded. Template
A's dietary runs were regenerated after the scratch allocation was reaped, two of twenty
peptides moved roughly ten times seed jitter, and the obvious explanation -- a rebuilt
environment -- could not be checked, because the only version recorded was align-it's.

These tests pin the two properties that matter: the record names the libraries that can
change a score, and a missing version degrades the record instead of failing a screen
that has already spent hours of CPU.
"""
import importlib.util
import sys
from pathlib import Path

import pytest

_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_ROOT / "src"))

_SPEC = importlib.util.spec_from_file_location(
    "screen_step", _ROOT / "scripts" / "04_screen.py")
screen = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(screen)


def test_record_names_every_library_that_can_change_a_score(monkeypatch):
    monkeypatch.setattr(screen, "version", lambda: "align-it v1.0.4")
    rec = screen.environment_record()
    # rdkit embeds the conformers, align-it scores them, obabel protonates between the
    # two. All three can move a number, so all three are recorded.
    for key in ("rdkit", "alignit", "obabel", "python", "platform"):
        assert key in rec, f"{key} missing from the provenance record"
    assert rec["rdkit"], "rdkit version must be a real value, not empty"


def test_a_missing_obabel_does_not_fail_the_screen(monkeypatch):
    """The screen has already spent its CPU by the time this is written."""
    def boom(*_a, **_k):
        raise FileNotFoundError("obabel")
    monkeypatch.setattr(screen.subprocess, "run", boom)
    monkeypatch.setattr(screen, "version", lambda: "align-it v1.0.4")
    rec = screen.environment_record()
    assert rec["obabel"] is None
    assert rec["rdkit"]


def test_an_unknown_version_is_recorded_as_null_not_omitted(monkeypatch):
    """A visible gap beats an absent key: absence is indistinguishable from an old run."""
    fake = type(sys)("rdkit")          # a module object with no __version__
    monkeypatch.setitem(sys.modules, "rdkit", fake)
    monkeypatch.setattr(screen, "version", lambda: "align-it v1.0.4")
    rec = screen.environment_record()
    assert "rdkit" in rec
    assert rec["rdkit"] is None


def test_the_record_is_json_serialisable():
    """It is written straight into block JSON; a non-serialisable value loses the run."""
    import json
    json.dumps(screen.environment_record())


@pytest.mark.parametrize("key", ["rdkit", "alignit", "obabel"])
def test_no_library_is_silently_dropped(key, monkeypatch):
    monkeypatch.setattr(screen, "version", lambda: "align-it v1.0.4")
    assert key in screen.environment_record()
