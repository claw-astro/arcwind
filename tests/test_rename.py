"""The old name ``diskwind`` still imports, as the arcwind modules themselves."""
import subprocess
import sys

import pytest


def test_diskwind_alias_is_arcwind():
    import arcwind
    import arcwind.fitting
    with pytest.warns(FutureWarning, match="renamed 'arcwind'"):
        import diskwind
    import diskwind.fitting
    from diskwind.fitting import FitResult
    assert diskwind is arcwind
    assert diskwind.fitting is arcwind.fitting
    assert FitResult is arcwind.FitResult


def test_legacy_command_runs():
    out = subprocess.run([sys.executable, "-c",
                          "from arcwind.cli import main_legacy; main_legacy(['options', 'fit'])"],
                         capture_output=True, text=True)
    assert out.returncode == 0
    assert "now 'arcwind'" in out.stderr
