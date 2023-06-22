# NOTE(hadim): typing_extensions can be replaced by typing once we drop support for Python 3.9.
from typing_extensions import TypeAlias
from typing import Union
from typing import Tuple

from rdkit import Chem
from rdkit.Chem import rdChemReactions

Mol: TypeAlias = Chem.rdchem.Mol
BondType: TypeAlias = Chem.rdchem.BondType
ChemicalReaction: TypeAlias = rdChemReactions.ChemicalReaction
Atom: TypeAlias = Chem.rdchem.Atom
Bond: TypeAlias = Chem.rdchem.Bond

RDKitColor = Union[Tuple[float, float, float, float], Tuple[float, float, float]]
DatamolColor = Union[RDKitColor, str]
