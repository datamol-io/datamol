from .errors import MolVSError as MolVSError, StandardizeError as StandardizeError, ValidateError as ValidateError
from .standardize import Standardizer as Standardizer, canonicalize_tautomer_smiles as canonicalize_tautomer_smiles, enumerate_tautomers_smiles as enumerate_tautomers_smiles, standardize_smiles as standardize_smiles
from .validate import Validator as Validator, validate_smiles as validate_smiles
from rdkit import Chem as Chem
from rdkit.Chem.MolStandardize import rdMolStandardize as rdMolStandardize
from typing import Any

log: Any

def ReorderTautomers(molecule: Any): ...
