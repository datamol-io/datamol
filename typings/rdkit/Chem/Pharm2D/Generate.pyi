from rdkit.Chem.Pharm2D import SigFactory as SigFactory, Utils as Utils
from rdkit.RDLogger import logger as logger
from typing import Any, Optional

def Gen2DFingerprint(mol: Any, sigFactory: Any, perms: Optional[Any] = ..., dMat: Optional[Any] = ..., bitInfo: Optional[Any] = ...): ...
