from rdkit import Chem as Chem
from rdkit.VLib.Filter import FilterNode as FilterNode
from typing import Any

class SmartsFilter(FilterNode):
    def __init__(self, patterns: Any = ..., counts: Any = ..., **kwargs: Any) -> None: ...
    def filter(self, cmpd: Any): ...
