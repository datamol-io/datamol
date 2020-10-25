from typing import Union
from typing import List
from typing import Tuple

from rdkit import Chem
from rdkit.Chem import Draw


def to_image(
    mols: Union[List[Chem.Mol], Chem.Mol],
    legends: Union[List[str], str] = None,
    n_cols: int = 4,
    use_svg: bool = False,
    mol_size: Tuple[int, int] = (300, 300),
    highlight_atom: List[List[int]] = None,
    highlight_bond: List[List[int]] = None,
):

    if isinstance(mols, Chem.Mol):
        mols = [mols]

    if isinstance(legends, str):
        legends = [legends]

    _highlight_atom = highlight_atom
    if highlight_atom is not None and isinstance(highlight_atom[0], int):
        _highlight_atom = [highlight_atom]

    _highlight_bond = highlight_bond
    if highlight_bond is not None and isinstance(highlight_bond[0], int):
        _highlight_bond = [highlight_bond]

    # Don't make the image bigger than it
    if len(mols) < n_cols:
        n_cols = len(mols)

    image = Draw.MolsToGridImage(
        mols,
        legends=legends,
        molsPerRow=n_cols,
        useSVG=use_svg,
        subImgSize=mol_size,
        highlightAtomLists=_highlight_atom,
        highlightBondLists=_highlight_bond,
    )
    return image
