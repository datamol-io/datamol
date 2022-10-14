__version__ = "0.7.17"


import rdkit
import packaging.version


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
