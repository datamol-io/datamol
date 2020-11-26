from rdkit import Chem as Chem
from rdkit.Chem import AllChem as AllChem, Crippen as Crippen, Descriptors as Descriptors, Lipinski as Lipinski
from rdkit.Dbase import DbModule as DbModule
from rdkit.Dbase.DbConnection import DbConnect as DbConnect
from typing import Any, Optional

logger: Any

def ProcessMol(mol: Any, typeConversions: Any, globalProps: Any, nDone: Any, nameProp: str = ..., nameCol: str = ..., redraw: bool = ..., keepHs: bool = ..., skipProps: bool = ..., addComputedProps: bool = ..., skipSmiles: bool = ..., uniqNames: Optional[Any] = ..., namesSeen: Optional[Any] = ...): ...
def ConvertRows(rows: Any, globalProps: Any, defaultVal: Any, skipSmiles: Any) -> None: ...
def LoadDb(suppl: Any, dbName: Any, nameProp: str = ..., nameCol: str = ..., silent: bool = ..., redraw: bool = ..., errorsTo: Optional[Any] = ..., keepHs: bool = ..., defaultVal: str = ..., skipProps: bool = ..., regName: str = ..., skipSmiles: bool = ..., maxRowsCached: int = ..., uniqNames: bool = ..., addComputedProps: bool = ..., lazySupplier: bool = ..., startAnew: bool = ...) -> None: ...
