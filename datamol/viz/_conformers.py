from typing import Union
from typing import List

import copy
import itertools

from rdkit import Chem
from rdkit.Chem import rdMolAlign


def _get_nglview():
    try:
        import nglview as nv

        return nv
    except ImportError:
        raise ImportError("You must install nglview from https://github.com/nglviewer/nglview.")


def _get_ipywidgets():
    try:
        import ipywidgets as widgets

        return widgets
    except ImportError:
        raise ImportError(
            "You must install ipywidgets from https://github.com/jupyter-widgets/ipywidgets/."
        )


def conformers(
    mol: Chem.rdchem.Mol,
    conf_id: int = -1,
    n_confs: Union[int, List[int]] = None,
    align_conf: bool = True,
    n_cols: int = 3,
    sync_views: bool = True,
    remove_hs: bool = True,
    width: str = "auto",
):
    """Visualize the conformer(s) of a molecule.

    Args:
        mol: a molecule.
        conf_id: The ID of the conformer to show. -1 shows
            the first conformer. Only works if `n_confs` is None.
        n_confs: Can be a number of conformers
            to shows or a list of conformer indices. When None, only the first
            conformer is displayed. When -1, show all conformers.
        align_conf: Whether to align conformers together.
        n_cols: Number of columns. Defaults to 3.
        sync_views: Wether to sync the multiple views.
        remove_hs: Wether to remove the hydrogens of the conformers.
        width: The width of the returned view. Defaults to "auto".
    """

    widgets = _get_ipywidgets()
    nv = _get_nglview()

    if mol.GetNumConformers() == 0:
        raise ValueError(
            "The molecule has 0 conformers. You can generate conformers with `dm.conformers.generate(mol)`."
        )

    # Clone the molecule
    mol = copy.deepcopy(mol)

    if remove_hs:
        mol = Chem.RemoveHs(mol)  # type: ignore
    else:
        mol = Chem.AddHs(mol)  # type: ignore

    if n_confs is None:
        return nv.show_rdkit(mol, conf_id=conf_id)

    # If n_confs is int, convert to list of conformer IDs
    if n_confs == -1:
        n_confs = [conf.GetId() for conf in mol.GetConformers()]
    elif isinstance(n_confs, int):
        if n_confs > mol.GetNumConformers():
            n_confs = mol.GetNumConformers()
        n_confs = list(range(n_confs))  # type: ignore

    if align_conf:
        rdMolAlign.AlignMolConformers(mol, confIds=n_confs)

    # Get number of rows
    n_rows = len(n_confs) // n_cols
    n_rows += 1 if (len(n_confs) % n_cols) > 0 else 0

    # Create a grid
    grid = widgets.GridspecLayout(n_rows, n_cols)  # type: ignore

    # Create and add views to the grid.
    widget_coords = itertools.product(range(n_rows), range(n_cols))
    views = []
    for i, (conf_id, (x, y)) in enumerate(zip(n_confs, widget_coords)):
        view = nv.show_rdkit(mol, conf_id=conf_id)
        view.layout.width = width
        view.layout.align_self = "stretch"
        grid[x, y] = view
        views.append(view)

    # Sync views
    if sync_views:
        for view in views:
            view._set_sync_camera(views)

    return grid
