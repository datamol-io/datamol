from rdkit import Chem as Chem, DataStructs as DataStructs
from rdkit.Chem import rdqueries as rdqueries
from typing import Any

rdkitFpParams: Any
FTYPE_ACYCLIC: str
FTYPE_CYCLIC: str
FTYPE_CYCLIC_ACYCLIC: str
ACYC_SMARTS: Any
CYC_SMARTS: Any
cSma1: Any
cSma2: Any
dummyAtomQuery: Any

def delete_bonds(mol: Any, bonds: Any, ftype: Any, hac: Any): ...
def select_fragments(fragments: Any, ftype: Any, hac: Any): ...
def isValidRingCut(mol: Any): ...
def generate_fraggle_fragmentation(mol: Any, verbose: bool = ...): ...
def atomContrib(subs: Any, mol: Any, tverskyThresh: float = ...): ...

modified_query_fps: Any

def compute_fraggle_similarity_for_subs(inMol: Any, qMol: Any, qSmi: Any, qSubs: Any, tverskyThresh: float = ...): ...
def GetFraggleSimilarity(queryMol: Any, refMol: Any, tverskyThresh: float = ...): ...
