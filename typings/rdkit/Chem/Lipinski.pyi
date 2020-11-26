from rdkit import Chem as Chem
from rdkit.Chem import rdMolDescriptors as rdMolDescriptors
from typing import Any

HDonorSmarts: Any
HAcceptorSmarts: Any
HeteroatomSmarts: Any
RotatableBondSmarts: Any
NHOHSmarts: Any
NOCountSmarts: Any
NumHDonors: Any
NumHAcceptors: Any
NumHeteroatoms: Any
NumRotatableBonds: Any
NOCount: Any
NHOHCount: Any
RingCount: Any

def HeavyAtomCount(mol: Any): ...

nm: Any
