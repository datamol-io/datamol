from rdkit import Geometry as Geometry, RDConfig as RDConfig
from typing import Any, Optional

logger: Any
BEGIN: int
END: int
TRIANGLE_FAN: int
COLOR: int
VERTEX: int
NORMAL: int
SPHERE: int
CYLINDER: int
ALPHA: int

def ShowArrow(viewer: Any, tail: Any, head: Any, radius: Any, color: Any, label: Any, transparency: int = ..., includeArrowhead: bool = ...) -> None: ...
def ShowMolFeats(mol: Any, factory: Any, viewer: Any, radius: float = ..., confId: int = ..., showOnly: bool = ..., name: str = ..., transparency: float = ..., colors: Optional[Any] = ..., excludeTypes: Any = ..., useFeatDirs: bool = ..., featLabel: Optional[Any] = ..., dirLabel: Optional[Any] = ..., includeArrowheads: bool = ..., writeFeats: bool = ..., showMol: bool = ..., featMapFile: bool = ...) -> None: ...

parser: Any
