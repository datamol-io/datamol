from rdkit import Chem as Chem
from rdkit.Chem import AllChem as AllChem
from typing import Any, Optional

def GetHeterocycleReactionSmarts(): ...

REACTION_CACHE: Any

def GetHeterocycleReactions(): ...
def EnumerateHeterocycles(inputmol: Any, depth: Optional[Any] = ...) -> None: ...
