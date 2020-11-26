from rdkit import Chem as Chem
from rdkit.Chem import Crippen as Crippen, rdMolDescriptors as rdMolDescriptors, rdPartialCharges as rdPartialCharges
from typing import Any, Optional

ptable: Any
bondScaleFacts: Any
mrBins: Any

def pySMR_VSA_(mol: Any, bins: Optional[Any] = ..., force: int = ...): ...

SMR_VSA_: Any
logpBins: Any

def pySlogP_VSA_(mol: Any, bins: Optional[Any] = ..., force: int = ...): ...

SlogP_VSA_: Any
chgBins: Any

def pyPEOE_VSA_(mol: Any, bins: Optional[Any] = ..., force: int = ...): ...

PEOE_VSA_: Any

def pyLabuteASA(mol: Any, includeHs: int = ...): ...

LabuteASA: Any
TPSA: Any
