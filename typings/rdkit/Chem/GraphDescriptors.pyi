from rdkit import Chem as Chem
from rdkit.Chem import Graphs as Graphs, rdMolDescriptors as rdMolDescriptors, rdchem as rdchem
from rdkit.ML.InfoTheory import entropy as entropy
from typing import Any, Optional

ptable: Any
hallKierAlphas: Any

def Ipc(mol: Any, avg: int = ..., dMat: Optional[Any] = ..., forceDMat: int = ...): ...

HallKierAlpha: Any
Kappa1: Any
Kappa2: Any
Kappa3: Any

def Chi0(mol: Any): ...
def Chi1(mol: Any): ...

Chi0v: Any
Chi1v: Any
Chi2v: Any
Chi3v: Any
Chi4v: Any
ChiNv_: Any
Chi0n: Any
Chi1n: Any
Chi2n: Any
Chi3n: Any
Chi4n: Any
ChiNn_: Any

def BalabanJ(mol: Any, dMat: Optional[Any] = ..., forceDMat: int = ...): ...
def BertzCT(mol: Any, cutoff: int = ..., dMat: Optional[Any] = ..., forceDMat: int = ...): ...
