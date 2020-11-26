from rdkit import DataStructs as DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols as FingerprintMols, MolSimilarity as MolSimilarity
from rdkit.ML.Cluster import Murtagh as Murtagh
from typing import Any

message = FingerprintMols.message
error = FingerprintMols.error

def GetDistanceMatrix(data: Any, metric: Any, isSimilarity: int = ...): ...
def ClusterPoints(data: Any, metric: Any, algorithmId: Any, haveLabels: bool = ..., haveActs: bool = ..., returnDistances: bool = ...): ...
def ClusterFromDetails(details: Any): ...
