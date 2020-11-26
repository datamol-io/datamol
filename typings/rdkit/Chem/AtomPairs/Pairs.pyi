from rdkit import DataStructs as DataStructs
from rdkit.Chem import rdMolDescriptors as rdMolDescriptors
from rdkit.Chem.AtomPairs import Utils as Utils
from rdkit.Chem.rdMolDescriptors import GetAtomPairFingerprint as GetAtomPairFingerprint, GetHashedAtomPairFingerprint as GetHashedAtomPairFingerprint
from typing import Any, Optional

GetAtomPairFingerprintAsIntVect: Any
numPathBits: Any
numFpBits: Any
fpLen: Any

def pyScorePair(at1: Any, at2: Any, dist: Any, atomCodes: Optional[Any] = ..., includeChirality: bool = ...): ...
def ExplainPairScore(score: Any, includeChirality: bool = ...): ...
def GetAtomPairFingerprintAsBitVect(mol: Any): ...
