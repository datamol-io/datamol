from rdkit import Chem as Chem
from rdkit.VLib.Transform import TransformNode as TransformNode
from typing import Any

class SmartsRemover(TransformNode):
    def __init__(self, patterns: Any = ..., wholeFragments: int = ..., **kwargs: Any) -> None: ...
    def transform(self, cmpd: Any): ...

biggerTest: str
__test__: Any
