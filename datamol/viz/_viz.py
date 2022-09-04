from typing import Union
from typing import List
from typing import Tuple
from typing import Optional
from typing import Any

import fsspec

from rdkit.Chem import Draw

import PIL

import datamol as dm


from .utils import prepare_mol_for_drawing


def to_image(
    mols: Union[List[dm.Mol], dm.Mol],
    legends: Union[List[Union[str, None]], str, None] = None,
    n_cols: int = 4,
    use_svg: bool = True,
    mol_size: Union[Tuple[int, int], int] = (300, 300),
    highlight_atom: Optional[List[List[int]]] = None,
    highlight_bond: Optional[List[List[int]]] = None,
    outfile: Optional[str] = None,
    max_mols: int = 32,
    copy: bool = True,
    indices: bool = False,
    bond_indices: bool = False,
    bond_line_width: int = 2,
    stereo_annotations: bool = True,
    legend_fontsize: int = 16,
    kekulize: bool = True,
    align: Union[dm.Mol, str, bool] = False,
    **kwargs: Any,
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
        bond_line_width: The width of the bond lines.
        legend_fontsize: Font size for the legend.
        kekulize: Run kekulization routine on molecules. Skipped if fails.
        align: Whether to align the 2D coordinates of the molecules.
            - If set to True, align all molecules with `dm.align.auto_align_many()`.
            - If set to a molecule, it is used as a template for alignment with `dm.align.template_align()`.
            - If set to False, no alignment is performed.
            For a more custom alignment, we suggest using directly the module `dm.align` instead.
        **kwargs: Additional arguments to pass to the drawing function. See RDKit
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

    # Whether to align the molecules
    if isinstance(align, (dm.Mol, str)):
        mols = [dm.align.template_align(mol, template=align) for mol in mols]
    elif align is True:
        mols = dm.align.auto_align_many(mols)

    # Prepare molecules before drawing
    mols = [prepare_mol_for_drawing(mol, kekulize=kekulize) for mol in mols]

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
    draw_options.bondLineWidth = bond_line_width

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
                    f.write(image.encode())  # type: ignore
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
