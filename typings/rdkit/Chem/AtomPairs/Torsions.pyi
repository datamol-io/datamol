from rdkit.Chem import rdMolDescriptors as rdMolDescriptors
from rdkit.Chem.AtomPairs import Utils as Utils
from rdkit.Chem.rdMolDescriptors import GetHashedTopologicalTorsionFingerprint as GetHashedTopologicalTorsionFingerprint, GetTopologicalTorsionFingerprint as GetTopologicalTorsionFingerprint
from typing import Any, Optional

GetTopologicalTorsionFingerprintAsIntVect: Any

def pyScorePath(mol: Any, path: Any, size: Any, atomCodes: Optional[Any] = ...): ...
def ExplainPathScore(score: Any, size: int = ...): ...
def GetTopologicalTorsionFingerprintAsIds(mol: Any, targetSize: int = ...): ...
