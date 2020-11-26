from rdkit import Chem as Chem
from rdkit.Chem import AllChem as AllChem, Crippen as Crippen
from rdkit.Chem.ChemUtils.AlignDepict import AlignDepict as AlignDepict
from typing import Any, Optional

logger: Any

def Usage() -> None: ...

nDumped: int

def Explode(template: Any, sidechains: Any, outF: Any, autoNames: bool = ..., do3D: bool = ..., useTethers: bool = ...) -> None: ...
def MoveDummyNeighborsToBeginning(mol: Any, useAll: bool = ...): ...
def ConstructSidechains(suppl: Any, sma: Optional[Any] = ..., replace: bool = ..., useAll: bool = ...): ...
