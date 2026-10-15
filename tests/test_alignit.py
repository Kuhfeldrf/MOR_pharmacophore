import pytest

from mor_pharmacophore.alignit import alignit_help, alignit_path

pytestmark = pytest.mark.skipif(alignit_path() is None, reason="align-it not installed")


def test_alignit_runs():
    assert alignit_help().returncode == 0
