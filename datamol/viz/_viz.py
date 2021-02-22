import fsspec

from typing import Union
from typing import List
from typing import Tuple

from rdkit import Chem
from rdkit.Chem import Draw

import PIL.PngImagePlugin


def to_image(
    mols: Union[List[Chem.Mol], Chem.Mol],
    legends: Union[List[str], str] = None,
    n_cols: int = 4,
    use_svg: bool = False,
    mol_size: Tuple[int, int] = (200, 200),
    highlight_atom: List[List[int]] = None,
    highlight_bond: List[List[int]] = None,
    outfile: str = None,
    max_mols: int = 50,
):
    """Generate an image out of a molecule or a list of molecule."""

    if isinstance(mols, Chem.Mol):
        mols = [mols]

    if isinstance(legends, str):
        legends = [legends]

    if max_mols is not None:
        mols = mols[:max_mols]

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

    if outfile is not None:
        with fsspec.open(outfile, "wb") as f:
            if use_svg:
                if isinstance(image, str):
                    # in a terminal process
                    f.write(image.encode())
                else:
                    # in a jupyter kernel process
                    f.write(image.data.encode())
            else:
                if isinstance(image, PIL.PngImagePlugin.PngImageFile):
                    # in a terminal process
                    image.save(f)
                else:
                    # in a jupyter kernel process
                    f.write(image.data)

    return image
