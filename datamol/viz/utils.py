from typing import List
from typing import Optional

from rdkit.Chem import AllChem
from rdkit import Geometry
from rdkit.Chem import Draw

import datamol as dm


def align_2d_coordinates(mols: List[dm.Mol], pattern: dm.Mol = None, **mcs_args):
    """Align the cooridnates of a list of molecules given a pattern or their
    maximum common substructure if no pattern is provided.

    Args:
        mol: List of molecules.
        pattern: Pattern to align the molecules to. If none is provided, the
            MCS will be computed.
        mcs_args: Arguments for the MCS (only when `pattern` is set).
    """

    if pattern is None:
        mcs = dm.find_mcs_with_details(mols=mols, **mcs_args)

        if mcs.smartsString is None or mcs.smartsString == "":
            # Do nothing
            return

        pattern = dm.from_smarts(mcs.smartsString)

    # Match the pattern to the molecules
    matches = [mol.GetSubstructMatch(pattern) for mol in mols]

    # Use the first molecule as a reference
    ref_mol, ref_match = mols[0], matches[0]
    AllChem.Compute2DCoords(ref_mol)  # type: ignore

    # Compute the coordinates of the ref molecule for the matching atoms
    coords = [ref_mol.GetConformer().GetAtomPosition(x) for x in ref_match]
    coords = [Geometry.Point2D(pt.x, pt.y) for pt in coords]  # type: ignore

    # Compute coords for the other molecules using
    # the pre-computed coords of the reference molecule.
    for mol, match in zip(mols[1:], matches[1:]):
        if not match:
            continue
        coord_dict = {match[i]: coord for i, coord in enumerate(coords)}
        AllChem.Compute2DCoords(mol, coordMap=coord_dict)  # type: ignore


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
