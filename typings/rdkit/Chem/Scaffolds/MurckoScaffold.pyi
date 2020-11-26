from rdkit import Chem as Chem
from rdkit.Chem import AllChem as AllChem
from typing import Any, Optional

murckoTransforms: Any

def MakeScaffoldGeneric(mol: Any): ...

murckoPatts: Any
murckoQ: Any
aromaticNTransform: Any

def GetScaffoldForMol(mol: Any): ...
def MurckoScaffoldSmiles(smiles: Optional[Any] = ..., mol: Optional[Any] = ..., includeChirality: bool = ...): ...
def MurckoScaffoldSmilesFromSmiles(smiles: Any, includeChirality: bool = ...): ...
