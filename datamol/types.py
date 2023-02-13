# NOTE(hadim): typing_extensions can be replaced by typing once we drop support for Python 3.9.
from typing_extensions import TypeAlias

from rdkit import Chem
from rdkit.Chem import rdChemReactions

Mol: TypeAlias = Chem.rdchem.Mol
BondType: TypeAlias = Chem.rdchem.BondType
ChemicalReaction: TypeAlias = rdChemReactions.ChemicalReaction
Atom: TypeAlias = Chem.rdchem.Atom
Bond: TypeAlias = Chem.rdchem.Bond
