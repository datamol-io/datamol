from rdkit.ML.Cluster import Clusters as Clusters
from rdkit.ML.Cluster.Clustering import MurtaghCluster as MurtaghCluster, MurtaghDistCluster as MurtaghDistCluster
from typing import Any

WARDS: int
SLINK: int
CLINK: int
UPGMA: int
MCQUITTY: int
GOWER: int
CENTROID: int
methods: Any

def ClusterData(data: Any, nPts: Any, method: Any, isDistData: int = ...): ...
