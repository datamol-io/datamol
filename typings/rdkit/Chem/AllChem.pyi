from rdkit.Chem import *
from rdkit.Chem.ChemicalFeatures import *
from rdkit.Chem.rdChemReactions import *
from rdkit.Chem.rdDepictor import *
from rdkit.Chem.rdDistGeom import *
from rdkit.Chem.rdForceFieldHelpers import *
from rdkit.Chem.rdMolAlign import *
from rdkit.Chem.rdMolDescriptors import *
from rdkit.Chem.rdMolTransforms import *
from rdkit.Chem.rdPartialCharges import *
from rdkit.Chem.rdReducedGraphs import *
from rdkit.Chem.rdShapeHelpers import *
from rdkit.Chem.rdqueries import *
from rdkit.Chem.rdMolEnumerator import *
from rdkit.Chem.rdSLNParse import *
from rdkit import DataStructs as DataStructs, ForceField as ForceField, RDConfig as RDConfig, rdBase as rdBase
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers as EnumerateStereoisomers, StereoEnumerationOptions as StereoEnumerationOptions
from rdkit.Geometry import rdGeometry as rdGeometry
from rdkit.RDLogger import logger as logger
from typing import Any, Optional

def TransformMol(mol: Any, tform: Any, confId: int = ..., keepConfs: bool = ...) -> None: ...
def ComputeMolShape(mol: Any, confId: int = ..., boxDim: Any = ..., spacing: float = ..., **kwargs: Any): ...
def ComputeMolVolume(mol: Any, confId: int = ..., gridSpacing: float = ..., boxMargin: float = ...): ...
def GetConformerRMS(mol: Any, confId1: Any, confId2: Any, atomIds: Optional[Any] = ..., prealigned: bool = ...): ...
def GetConformerRMSMatrix(mol: Any, atomIds: Optional[Any] = ..., prealigned: bool = ...): ...
def EnumerateLibraryFromReaction(reaction: Any, sidechainSets: Any, returnReactants: bool = ...) -> None: ...
def ConstrainedEmbed(mol: Any, core: Any, useTethers: bool = ..., coreConfId: int = ..., randomseed: int = ..., getForceField: Any = ..., **kwargs: Any): ...
def AssignBondOrdersFromTemplate(refmol: Any, mol: Any): ...
