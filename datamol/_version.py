try:
    from importlib.metadata import version
    from importlib.metadata import PackageNotFoundError
except ModuleNotFoundError:
    # Try backported to PY<38 `importlib_metadata`.
    from importlib_metadata import version
    from importlib_metadata import PackageNotFoundError


import rdkit
import packaging.version


try:
    __version__ = version("datamol")
except PackageNotFoundError:
    # package is not installed
    __version__ = "dev"

CURRENT_RDKIT_VERSION = rdkit.__version__
CURRENT_RDKIT_VERSION_OBJ = packaging.version.parse(CURRENT_RDKIT_VERSION)


def is_lower_than_current_rdkit_version(rdkit_version: str):
    return CURRENT_RDKIT_VERSION_OBJ < packaging.version.parse(rdkit_version)


def is_greater_than_current_rdkit_version(rdkit_version: str):
    return CURRENT_RDKIT_VERSION_OBJ > packaging.version.parse(rdkit_version)


def is_lower_eq_than_current_rdkit_version(rdkit_version: str):
    return CURRENT_RDKIT_VERSION_OBJ <= packaging.version.parse(rdkit_version)


def is_greater_eq_than_current_rdkit_version(rdkit_version: str):
    return CURRENT_RDKIT_VERSION_OBJ >= packaging.version.parse(rdkit_version)
