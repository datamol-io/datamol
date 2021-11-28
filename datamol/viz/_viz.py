import itertools

import fsspec

from typing import Union
from typing import List
from typing import Tuple

from rdkit import Chem
from rdkit.Chem import Draw

import PIL

import datamol as dm


def to_image(
    mols: Union[List[dm.Mol], dm.Mol],
    legends: Union[List[Union[str, None]], str, None] = None,
    n_cols: int = 4,
    use_svg: bool = False,
    mol_size: Union[Tuple[int, int], int] = (200, 200),
    highlight_atom: List[List[int]] = None,
    highlight_bond: List[List[int]] = None,
    outfile: str = None,
    max_mols: int = 32,
    copy: bool = True,
    indices: bool = False,
    bond_indices: bool = False,
    stereo_annotations: bool = True,
    legend_fontsize: int = 16,
    kekulize: bool = True,
    align: Union[bool, dm.Mol] = False,
    **kwargs,
):
    """Generate an image out of a molecule or a list of molecules.

    Args:
        mols: One or a list of molecules.
        legends: A string or a list of string as legend for every molecules.
        n_cols: Number of molecules per column.
        use_svg: Whether to ouput an SVG (or a PNG).
        mol_size: A int or a tuple of int defining the size per molecule.
        highlight_atom: the atoms to highlight.
        highlight_bond: The bonds to highlight.
        outfile: Path where to save the image (local or remote path).
        max_mols: The maximum number of molecules to display.
        copy: Whether to copy the molecules or not.
        indices: Whether to draw the atom indices.
        bond_indices: Whether to draw the bond indices.
        legend_fontsize: Font size for the legend.
        kekulize: Run kekulization routine on molecules. Skipped if fails.
        align: Whether to align the 2D coordinates of the molecules. If True
            or set to a valid molecule object `dm.viz.utils.align_2d_coordinates` is used.
            If `align` is set to a molecule object, this molecule will be used as a pattern
            for the alignment. If `align` is set to True, the MCS will be computed.
        kwargs: Additional arguments to pass to the drawing function. See RDKit
            documentation related to `MolDrawOptions` for more details at
            https://www.rdkit.org/docs/source/rdkit.Chem.Draw.rdMolDraw2D.html.
    """

    if isinstance(mol_size, int):
        mol_size = (mol_size, mol_size)

    if isinstance(mols, dm.Mol):
        mols = [mols]

    if isinstance(legends, str):
        legends = [legends]

    if copy:
        mols = [dm.copy_mol(mol) for mol in mols]

    if max_mols is not None:
        mols = mols[:max_mols]

        if legends is not None:
            legends = legends[:max_mols]

    # Prepare molecules before drawing
    # Code is inspired from `rdkit.Chem.Draw._moltoimg`.
    _mols = []
    for mol in mols:
        if mol is not None:
            try:
                with dm.without_rdkit_log():
                    try:
                        mol.GetAtomWithIdx(0).GetExplicitValence()  # type: ignore
                    except RuntimeError:
                        mol.UpdatePropertyCache(False)  # type: ignore
                    _kekulize = Draw._okToKekulizeMol(mol, kekulize)
                    _mol = Draw.rdMolDraw2D.PrepareMolForDrawing(mol, kekulize=_kekulize)
            except ValueError:  # <- can happen on a kekulization failure
                _mol = Draw.rdMolDraw2D.PrepareMolForDrawing(mol, kekulize=False)
        else:
            _mol = None

        _mols.append(_mol)
    mols = _mols

    # Whether to align the molecules
    if align is True:
        dm.viz.utils.align_2d_coordinates(mols)
    elif isinstance(align, dm.Mol):
        dm.viz.utils.align_2d_coordinates(mols, pattern=align)

    _highlight_atom = highlight_atom
    if highlight_atom is not None and isinstance(highlight_atom[0], int):
        _highlight_atom = [highlight_atom]

    _highlight_bond = highlight_bond
    if highlight_bond is not None and isinstance(highlight_bond[0], int):
        _highlight_bond = [highlight_bond]

    # Don't make the image bigger than it
    if len(mols) < n_cols:
        n_cols = len(mols)

    draw_options = Draw.rdMolDraw2D.MolDrawOptions()
    draw_options.legendFontSize = legend_fontsize
    draw_options.addAtomIndices = indices
    draw_options.addBondIndices = bond_indices
    draw_options.addStereoAnnotation = stereo_annotations

    # Add the custom drawing options.
    _kwargs = {}
    for k, v in kwargs.items():
        if hasattr(draw_options, k):
            setattr(draw_options, k, v)
        else:
            _kwargs[k] = v

    image = Draw.MolsToGridImage(
        mols,
        legends=legends,
        molsPerRow=n_cols,
        useSVG=use_svg,
        subImgSize=mol_size,
        highlightAtomLists=_highlight_atom,
        highlightBondLists=_highlight_bond,
        drawOptions=draw_options,
        **_kwargs,
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
