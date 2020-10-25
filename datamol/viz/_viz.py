import copy
import itertools

import ipywidgets as widgets

from rdkit import Chem


def _get_nglview():
    try:
        import nglview as nv

        return nv
    except ImportError:
        raise ImportError("You must install nglview from https://github.com/nglviewer/nglview.")


def conformers(
    mol: Chem.Mol,
    conf_id: int = -1,
    n_confs=None,
    align_conf: bool = True,
    n_cols: int = 3,
    sync_views: bool = True,
    remove_hs: bool = True,
    width: str = "auto",
):
    """Visualize conformer(s) of a molecule.

    Args:
        mol (Chem.Mol): [description]
        conf_id (int, optional): [description]. Defaults to -1.
        n_confs ([type], optional): [description]. Defaults to None.
        align_conf (bool, optional): [description]. Defaults to True.
        n_cols (int, optional): [description]. Defaults to 3.
        sync_views (bool, optional): [description]. Defaults to True.
        remove_hs (bool, optional): [description]. Defaults to True.
        width (str, optional): [description]. Defaults to "50%".

    Returns:
        [type]: [description]
    """

    nv = _get_nglview()

    # Clone the molecule
    mol = copy.deepcopy(mol)

    if remove_hs:
        mol = Chem.RemoveHs(mol)
    else:
        mol = Chem.AddHs(mol)

    if n_confs is None:
        return nv.show_rdkit(mol, conf_id=conf_id)

    # If n_confs is int, convert to list of conformer IDs
    if isinstance(n_confs, int):
        if n_confs > mol.GetNumConformers():
            n_confs = mol.GetNumConformers()
        n_confs = list(range(n_confs))

    if align_conf:
        Chem.rdMolAlign.AlignMolConformers(mol, confIds=n_confs)

    # Get number of rows
    n_rows = len(n_confs) // n_cols
    n_rows += 1 if (len(n_confs) % n_cols) > 0 else 0

    # Create a grid
    grid = widgets.GridspecLayout(n_rows, n_cols)

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
