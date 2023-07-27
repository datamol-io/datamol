from typing import Optional
from typing import Union

import io
import fsspec

from rdkit.Chem import Draw
from matplotlib import colors as mcolors

import PIL.Image
import PIL.PngImagePlugin

import datamol as dm

from datamol.types import RDKitColor
from datamol.types import DatamolColor


def prepare_mol_for_drawing(mol: Optional[dm.Mol], kekulize: bool = True) -> Optional[dm.Mol]:
    """Prepare the molecule before drawing to avoid any error due to unsanitized molecule
    or incorrect valence or aromaticity.

    Code is inspired from `rdkit.Chem.Draw._moltoimg`.

    Args:
        mol: A molecule to prepare. If set to None, the function will return None.
        kekulize: Whether to kekulize the molecule.
    """

    if mol is None:
        return None

    try:
        with dm.without_rdkit_log():
            # Check for implicit and explicit valence
            if mol.NeedsUpdatePropertyCache():  # type: ignore
                mol.UpdatePropertyCache(False)  # type: ignore

            # Check for aromaticity
            if dm.is_lower_than_current_rdkit_version("2022.09"):
                _kekulize = Draw._okToKekulizeMol(mol, kekulize)  # type: ignore
            else:
                _kekulize = Draw.shouldKekulize(mol, kekulize)

            # Run the rdkit preparation procedure
            _mol = Draw.rdMolDraw2D.PrepareMolForDrawing(mol, kekulize=_kekulize)

    except ValueError:  # <- can happen on a kekulization failure
        # Run the rdkit preparation procedure with kekulize set to `False`
        _mol = Draw.rdMolDraw2D.PrepareMolForDrawing(mol, kekulize=False)

    return _mol


def is_ipython_session() -> bool:
    try:
        kernel_name = get_ipython().__class__.__name__  # noqa: F821 # type: ignore
        module_name = get_ipython().__class__.__module__  # noqa: F821 # type: ignore

        if kernel_name == "ZMQInteractiveShell" or module_name == "google.colab._shell":
            return True
    except Exception:
        pass

    return False


def drawer_to_image(drawer: Draw.rdMolDraw2D.MolDraw2D):
    """Convert an RDkit drawer to an image. The image can be either a PNG or SVG depending on the
    drawer class. The returned image type will depends whether the Python session is an IPython one or not.

    This function matches the behavior of `datamol.to_image` and `rdkit.Chem.Draw.MolDraw2DToImage`.

    Args:
        drawer: An RDkit drawer.

    Returns:
        An image: either PNG or SVG depending on the drawer class. If within an IPython sessions,
                  IPython display objects are returned.
    """

    is_svg = isinstance(drawer, Draw.rdMolDraw2D.MolDraw2DSVG)

    if is_ipython_session():
        if is_svg:
            from IPython.core.display import SVG

            return SVG(drawer.GetDrawingText())
        else:
            from IPython.core.display import Image

            return Image(drawer.GetDrawingText())
    else:
        if is_svg:
            return drawer.GetDrawingText()
        else:
            from PIL import Image

            return Image.open(io.BytesIO(drawer.GetDrawingText()))


def image_to_file(
    image: Union[
        str,
        PIL.PngImagePlugin.PngImageFile,
        bytes,
        PIL.Image.Image,
    ],
    outfile,
    as_svg: bool = False,
):
    """Save image to file. The image can be either a PNG or SVG depending

    Args:
        image: Image to save to a file
        outfile: Path to the output file where to save the image
        as_svg: Whether the image is an SVG or not
    """

    with fsspec.open(outfile, "wb") as f:
        if as_svg:
            if isinstance(image, str):
                # in a terminal process
                f.write(image.encode())  # type: ignore
            else:
                # in a jupyter kernel process
                f.write(image.data.encode())  # type: ignore
        else:
            if isinstance(image, PIL.PngImagePlugin.PngImageFile):  # type: ignore
                # in a terminal process
                image.save(f)  # type: ignore
            else:
                # in a jupyter kernel process
                f.write(image.data)  # type: ignore


def to_rdkit_color(color: Optional[DatamolColor]) -> Optional[RDKitColor]:
    """If required convert a datamol color (rgb, rgba or hex string) to an RDKit
    color (rgb or rgba).

    Args:
        color: A datamol color: hex, rgb, rgba or None.
    """
    if isinstance(color, str):
        return mcolors.to_rgba(color)  # type: ignore
    return color
