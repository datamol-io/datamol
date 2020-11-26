from rdkit import Chem as Chem
from typing import Any

aldehydePatt: Any
ketonePatt: Any
amidePatt: Any
esterPatt: Any
carboxylatePatt: Any
carboxylPatt: Any
specialCases: Any

def SATISTypes(mol: Any, neighborsToInclude: int = ...): ...
