from typing import TypeAlias
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


def _get_explicit_valence(atom: Atom) -> int:
    """Return an atom's explicit valence across supported RDKit versions."""
    if hasattr(atom, "GetValence"):
        return atom.GetValence(Chem.ValenceType.EXPLICIT)
    return atom.GetExplicitValence()


def _get_implicit_valence(atom: Atom) -> int:
    """Return an atom's implicit valence across supported RDKit versions."""
    if hasattr(atom, "GetValence"):
        return atom.GetValence(Chem.ValenceType.IMPLICIT)
    return atom.GetImplicitValence()
