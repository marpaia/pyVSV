"""Check a built distribution in a fresh environment with Python's -I flag."""

import subprocess
import sys
import sysconfig
from importlib import import_module
from importlib.metadata import distribution
from importlib.resources import files
from pathlib import Path

import vsv


def main() -> None:
    package = distribution("pyVSV")
    assert package.metadata["Name"] == "pyVSV"
    assert package.version == vsv.__version__
    assert package.metadata["License-Expression"] == "Apache-2.0"
    assert package.metadata["Requires-Python"] == ">=3.12"
    assert any(str(path).endswith("licenses/LICENSE") for path in package.files or ())
    assert vsv.__file__ is not None
    assert Path(vsv.__file__).resolve().is_relative_to(Path(sys.prefix).resolve())
    assert not Path(vsv.__file__).resolve().is_relative_to(Path(__file__).resolve().parents[1])

    for module in ("vsv.cli", "vsv.design.sbol", "vsv.genome.sbol"):
        import_module(module)
    for resource in ("py.typed", "synthetic-virology.ttl", "NC_001560.ttl"):
        assert files("vsv").joinpath(resource).is_file(), f"Missing resource: {resource}"

    command = Path(sysconfig.get_path("scripts")) / "vsv"
    for entrypoint in ([sys.executable, "-I", "-m", "vsv"], [str(command)]):
        result = subprocess.run(
            [*entrypoint, "--version"], check=True, capture_output=True, text=True
        )
        assert result.stdout.strip() == f"vsv {package.version}"
        result = subprocess.run([*entrypoint, "--help"], check=True, capture_output=True, text=True)
        assert "Usage:" in result.stdout
    print(f"pyVSV {package.version}: installed package check passed")


if __name__ == "__main__":
    main()
