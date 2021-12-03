from typing import List
from typing import Optional

from rdkit.Chem import AllChem
from rdkit.Chem import rdDepictor
from rdkit.Chem import Draw

import datamol as dm


def align_2d_coordinates(
    mols: List[dm.Mol],
    pattern: dm.Mol = None,
    copy: bool = True,
    **mcs_args,
):
    """Align the cooridnates of a list of molecules given a pattern or their
    maximum common substructure if no pattern is provided.

    Args:
        mol: List of molecules.
        pattern: Pattern to align the molecules to. If none is provided, the
            MCS will be computed.
        copy: Whether to copy the molecules before aligning them.
        mcs_args: Arguments for the MCS (only when `pattern` is set).
    """

    if copy is True:
        mols = [dm.copy_mol(mol) for mol in mols]

    if pattern is None:
        mcs_smarts = dm.find_mcs(mols=mols, **mcs_args)

        if mcs_smarts is None:
            # Do nothing
            return mols

        pattern = dm.from_smarts(mcs_smarts)

    AllChem.Compute2DCoords(pattern)  # type: ignore

    for mol in mols:
        rdDepictor.GenerateDepictionMatching2DStructure(
            mol,
            reference=pattern,
            acceptFailure=True,
            allowRGroups=True,
        )

    return mols


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
