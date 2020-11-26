from rdkit.DataStructs.VectCollection import VectCollection as VectCollection
from rdkit.ML import InfoTheory as InfoTheory
from rdkit.ML.DecTree import SigTree as SigTree
from rdkit.ML.FeatureSelect import CMIM as CMIM
from typing import Any, Optional

def BuildSigTree(examples: Any, nPossibleRes: Any, ensemble: Optional[Any] = ..., random: int = ..., metric: Any = ..., biasList: Any = ..., depth: int = ..., maxDepth: int = ..., useCMIM: int = ..., allowCollections: bool = ..., verbose: int = ..., **kwargs: Any): ...
def SigTreeBuilder(examples: Any, attrs: Any, nPossibleVals: Any, initialVar: Optional[Any] = ..., ensemble: Optional[Any] = ..., randomDescriptors: int = ..., **kwargs: Any): ...
