from rdkit import Chem as Chem
from rdkit.Chem.Pharm2D import Utils as Utils
from typing import Any, Optional

class MatchError(Exception): ...

def GetAtomsMatchingBit(sigFactory: Any, bitIdx: Any, mol: Any, dMat: Optional[Any] = ..., justOne: int = ..., matchingAtoms: Optional[Any] = ...): ...
