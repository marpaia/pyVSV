"""Exercise the installed package through both supported entry points."""

import shutil
import subprocess
import sys
from importlib.metadata import version
from pathlib import Path

import pytest


@pytest.fixture(params=["console", "module"])
def command(request: pytest.FixtureRequest) -> list[str]:
    if request.param == "module":
        return [sys.executable, "-I", "-m", "vsv"]
    executable = shutil.which("vsv")
    assert executable is not None, "The installed vsv console script must be on PATH."
    return [executable]


def test_version(command: list[str], tmp_path: Path) -> None:
    result = subprocess.run(
        [*command, "--version"],
        cwd=str(tmp_path),
        capture_output=True,
        text=True,
        check=True,
    )
    assert result.stdout.strip() == f"vsv {version('vsv')}"
    assert result.stderr == ""


def test_help(command: list[str]) -> None:
    result = subprocess.run([*command, "--help"], capture_output=True, text=True, check=True)
    assert "usage: vsv" in result.stdout
    assert "--version" in result.stdout
    assert result.stderr == ""
