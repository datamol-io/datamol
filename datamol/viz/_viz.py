import fsspec

from typing import Union
from typing import List
from typing import Tuple

from rdkit import Chem
from rdkit.Chem import Draw

import PIL

import datamol as dm


def to_image(
    mols: Union[List[Chem.rdchem.Mol], Chem.rdchem.Mol],
    legends: Union[List[Union[str, None]], str, None] = None,
    n_cols: int = 4,
    use_svg: bool = False,
    mol_size: Union[Tuple[int, int], int] = (200, 200),
    highlight_atom: List[List[int]] = None,
    highlight_bond: List[List[int]] = None,
    outfile: str = None,
    max_mols: int = 32,
    copy: bool = False,
    indices: bool = False,
):
    """Generate an image out of a molecule or a list of molecule.

    Args:
        mols: one or a list of molecules.
        legends: a string or a list of string as legend for every molecules.
        n_cols: number of molecules per column.
        use_svg: whether to ouput an SVG (or a PNG).
        mol_size: a int or a tuple of int defining the size per molecule.
        highlight_atom: atom to highlight.
        highlight_bond: bonds to highlight.
        outfile: path where to save the image (local or remote path).
        max_mols: the maximum number of molecules to display.
        copy: whether to copy the molecules or not.
        indices: Whether to draw the atom indices.
    """

    if isinstance(mol_size, int):
        mol_size = (mol_size, mol_size)

    if isinstance(mols, Chem.rdchem.Mol):
        mols = [mols]

    if isinstance(legends, str):
        legends = [legends]

    if copy:
        mols = [dm.copy_mol(mol) for mol in mols]

    if max_mols is not None:
        mols = mols[:max_mols]

        if legends is not None:
            legends = legends[:max_mols]

    if indices is True:
        [dm.atom_indices_to_mol(mol) for mol in mols]

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
                    f.write(image.data.encode())  # type: ignore
            else:
                if isinstance(image, PIL.PngImagePlugin.PngImageFile):  # type: ignore
                    # in a terminal process
                    image.save(f)
                else:
                    # in a jupyter kernel process
                    f.write(image.data)  # type: ignore

    return image
