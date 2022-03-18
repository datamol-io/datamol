from typing import Optional

from rdkit.Chem import Draw

import datamol as dm


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
            _kekulize = Draw._okToKekulizeMol(mol, kekulize)

            # Run the rdkit preparation procedure
            _mol = Draw.rdMolDraw2D.PrepareMolForDrawing(mol, kekulize=_kekulize)

    except ValueError:  # <- can happen on a kekulization failure

        # Run the rdkit preparation procedure with kekulize set to `False`
        _mol = Draw.rdMolDraw2D.PrepareMolForDrawing(mol, kekulize=False)

    return _mol
