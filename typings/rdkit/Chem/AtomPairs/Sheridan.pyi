from rdkit import Chem as Chem, RDConfig as RDConfig
from rdkit.Chem import rdMolDescriptors as rdMolDescriptors
from rdkit.Chem.rdMolDescriptors import GetAtomPairFingerprint as GetAtomPairFingerprint, GetTopologicalTorsionFingerprint as GetTopologicalTorsionFingerprint
from typing import Any, Optional

numPathBits: Any
numFpBits: Any
fpLen: Any

def AssignPattyTypes(mol: Any, defns: Optional[Any] = ...): ...

typMap: Any

def GetBPFingerprint(mol: Any, fpfn: Any = ...): ...
def GetBTFingerprint(mol: Any, fpfn: Any = ...): ...
