"""Regression tests for the OpenBabel plugin-path collision that killed job 242423.

Two OpenBabels can coexist in one prefix and they are not interchangeable: align_it is
linked against conda-forge's libopenbabel, while the pip `openbabel` wheel that PLIP's
dependency chain installs bundles its own library and plugins built with a different C++
ABI. Importing the wheel sets BABEL_LIBDIR/BABEL_DATADIR into os.environ, the align_it
subprocess inherits them, and reading any structure then fails with

    undefined symbol: _ZN9OpenBabel8OBFormat7DisplayERSsPKcS3_

The nasty part is that `align_it --help` succeeds anyway, because help loads no format
plugin. A smoke test passes and the template build still dies. These tests therefore
assert on the environment handed to the subprocess, not on whether the binary runs.
"""
import os

import pytest

from mor_pharmacophore import alignit


@pytest.fixture
def fake_prefix(tmp_path, monkeypatch):
    """A prefix laid out like a conda env with align_it and OpenBabel 3.2.1 installed."""
    (tmp_path / "bin").mkdir()
    exe = tmp_path / "bin" / "align_it"
    exe.write_text("#!/bin/sh\nexit 0\n")
    exe.chmod(0o755)
    (tmp_path / "lib" / "openbabel" / "3.2.1").mkdir(parents=True)
    (tmp_path / "share" / "openbabel" / "3.2.1").mkdir(parents=True)
    monkeypatch.setattr(alignit, "alignit_path", lambda: str(exe))
    return tmp_path


def test_babel_dirs_are_pinned_to_the_binarys_own_prefix(fake_prefix):
    env = alignit._babel_env()
    assert env["BABEL_LIBDIR"] == str(fake_prefix / "lib" / "openbabel" / "3.2.1")
    assert env["BABEL_DATADIR"] == str(fake_prefix / "share" / "openbabel" / "3.2.1")


def test_a_polluted_babel_libdir_is_overridden_not_inherited(fake_prefix, monkeypatch):
    """This is the actual bug: the wheel's value must not reach the subprocess."""
    wheel = "/env/lib/python3.11/site-packages/openbabel/lib/openbabel/3.2.1"
    monkeypatch.setenv("BABEL_LIBDIR", wheel)
    monkeypatch.setenv("BABEL_DATADIR", wheel.replace("/lib/", "/share/"))
    env = alignit._babel_env()
    assert env["BABEL_LIBDIR"] != wheel
    assert "site-packages" not in env["BABEL_LIBDIR"]
    assert "site-packages" not in env["BABEL_DATADIR"]


def test_a_polluted_value_is_dropped_when_the_prefix_has_no_plugins(tmp_path, monkeypatch):
    """Inheriting a wrong plugin dir is worse than inheriting none at all."""
    (tmp_path / "bin").mkdir()
    exe = tmp_path / "bin" / "align_it"
    exe.write_text("#!/bin/sh\nexit 0\n")
    exe.chmod(0o755)
    monkeypatch.setattr(alignit, "alignit_path", lambda: str(exe))
    monkeypatch.setenv("BABEL_LIBDIR", "/somewhere/wrong")
    env = alignit._babel_env()
    assert "BABEL_LIBDIR" not in env


def test_rest_of_the_environment_is_preserved(fake_prefix, monkeypatch):
    monkeypatch.setenv("MOR_CANARY", "keep-me")
    env = alignit._babel_env()
    assert env["MOR_CANARY"] == "keep-me"
    assert "PATH" in env


def test_no_binary_means_the_environment_is_left_alone(monkeypatch):
    monkeypatch.setattr(alignit, "alignit_path", lambda: None)
    assert alignit._babel_env() == dict(os.environ)
