from rdkit import Chem as Chem
from typing import Any, Optional

environs: Any
reactionDefs: Any
smartsGps: Any
g1: Any
g2: Any
bnd: Any
r1: Any
r2: Any
sma: Any
t: Any
environMatchers: Any
bondMatchers: Any
tmp: Any
e1: Any
e2: Any
patt: Any
reactions: Any
reverseReactions: Any
rs: Any
ps: Any
rxn: Any
labels: Any

def FindBRICSBonds(mol: Any, randomizeOrder: bool = ..., silent: bool = ...) -> None: ...
def BreakBRICSBonds(mol: Any, bonds: Optional[Any] = ..., sanitize: bool = ..., silent: bool = ...): ...
def BRICSDecompose(mol: Any, allNodes: Optional[Any] = ..., minFragmentSize: int = ..., onlyUseReactions: Optional[Any] = ..., silent: bool = ..., keepNonLeafNodes: bool = ..., singlePass: bool = ..., returnMols: bool = ...): ...

dummyPattern: Any

def BRICSBuild(fragments: Any, onlyCompleteMols: bool = ..., seeds: Optional[Any] = ..., uniquify: bool = ..., scrambleReagents: bool = ..., maxDepth: int = ...) -> None: ...
