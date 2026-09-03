"""Exercise the installed distribution without importing the checkout."""

import os
from importlib.metadata import distribution
from pathlib import Path

import datamol as package

installed = distribution("datamol")
location = Path(package.__file__).resolve()
assert location == Path(installed.locate_file("datamol/__init__.py")).resolve(), location
assert not location.is_relative_to(Path(__file__).resolve().parents[2]), location
assert package.__version__ == installed.version
if expected := os.environ.get("EXPECTED_VERSION"):
    assert installed.version == expected, (installed.version, expected)

assert package.to_smiles(package.to_mol("CCO")) == "CCO"
print(f"datamol {installed.version}: {location}")
