from collections import namedtuple
from rdkit import Chem as Chem, RDConfig as RDConfig
from rdkit.Avalon import pyAvalonTools as pyAvalonTools
from rdkit.Chem.MolKey import InchiInfo as InchiInfo
from typing import Any, Optional

class MolIdentifierException(Exception): ...
class BadMoleculeException(Exception): ...

MOL_KEY_VERSION: str
ERROR_DICT: Any
INCHI_COMPUTATION_ERROR: Any
RDKIT_CONVERSION_ERROR: Any
INCHI_READWRITE_ERROR: Any
NULL_MOL: Any
BAD_SET: Any
GET_STEREO_RE: Any
NULL_SMILES_RE: Any
PATTERN_NULL_MOL: str
CHIRAL_POS: int
T_NULL_MOL: Any
stereo_code_dict: Any

def initStruchk(configDir: Optional[Any] = ..., logFile: Optional[Any] = ...) -> None: ...
def CheckCTAB(ctab: Any, isSmiles: bool = ...): ...

InchiResult = namedtuple('InchiResult', ['error', 'inchi', 'fixed_ctab'])

def GetInchiForCTAB(ctab: Any): ...
def ErrorBitsToText(err: Any): ...

MolKeyResult = namedtuple('MolKeyResult', ['mol_key', 'error', 'inchi', 'fixed_ctab', 'stereo_code', 'stereo_comment'])

def GetKeyForCTAB(ctab: Any, stereo_info: Optional[Any] = ..., stereo_comment: Optional[Any] = ..., logger: Optional[Any] = ...): ...
