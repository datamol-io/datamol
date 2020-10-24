from typing import Union

from rdkit import Chem


PERIODIC_TABLE = Chem.rdchem.GetPeriodicTable()


def to_mol(
    mol: Union[str, Chem.Mol],
    add_hs: bool = False,
    explicit_only: bool = True,
    ordered: bool = False,
    sanitize: bool = True,
):
    """Convert an input molecule (smiles representation) into a `Chem.Mol`

    Raise `ValueError` if the input is neither a Chem.Mol nor a string.

    Args:
        mol (str or rdkit.Chem.Mol): SMILES of a molecule or a molecule.
        add_hs (bool, optional): Whether hydrogens should be added the molecule. Default to False.
        explicit_only (bool, optional): Whether to only add explicit hydrogen or both
            (implicit and explicit). when `add_hs` is set to True. Default to True.
        ordered (bool, optional): Whether the atom should be ordered. This option is
            important if you want to ensure that the features returned will always maintain
            a single atom order for the same molecule, regardless of its original SMILES representation.
            Default to False.
        sanitize (bool optional): Whether to apply rdkit sanitization when input is a SMILES. Default to True.

    Returns:
        mol (rdkit.Chem.Mol): the molecule if some conversion have been made. If the conversion fails
        None is returned so make sure that you handle this case on your own.
    """

    if not isinstance(mol, (str, Chem.Mol)):
        raise ValueError(f"Input should be a Chem.Mol or a string instead of '{type(mol)}'")

    if isinstance(mol, str):
        mol = Chem.MolFromSmiles(mol, sanitize=sanitize)

        if not sanitize and mol is not None:
            mol.UpdatePropertyCache(False)

    # Add hydrogens
    if mol is not None and add_hs:
        mol = Chem.AddHs(mol, explicitOnly=explicit_only)

    # Reorder atoms
    if mol and ordered:
        mol = reorder_atoms(mol)

    return mol


def reorder_atoms(mol: Chem.Mol):
    """Reorder the atoms in a mol. It ensures a single atom order for the same molecule,
    regardless of its original representation.

    Args:
        mol (Chem.Mol): a molecule.

    Returns:
        mol (Chem.Mol): a molecule.
    """
    if mol.GetNumAtoms() == 0:
        return mol

    new_order = Chem.rdmolfiles.CanonicalRankAtoms(mol, breakTies=True)
    new_order = sorted([(y, x) for x, y in enumerate(new_order)])
    return Chem.RenumberAtoms(mol, [y for (x, y) in new_order])


def to_neutral(mol: Chem.Mol):
    """Neutralize the charge of a molecule.

    Args:
        mol (Chem.Mol): a molecule.

    Returns:
        mol (Chem.Mol): a molecule.
    """
    if mol is None:
        return mol

    for a in mol.GetAtoms():
        if a.GetFormalCharge() < 0 or (
            a.GetExplicitValence() >= PERIODIC_TABLE.GetDefaultValence(a.GetSymbol())
            and a.GetFormalCharge() > 0
        ):
            a.SetFormalCharge(0)
            a.UpdatePropertyCache(False)
    return mol


def sanitize_mol(mol: Chem.Mol, charge_neutral=False):
    """Sanitize molecule and fix common errors.

    Args:
        mol (Chem.Mol): a molecule.

    Returns:
        mol (Chem.Mol): a molecule.
        charge_neutral (bool, optional): whether charge neutralization should be applied

    Returns:
        mol or None
    """
    if charge_neutral:
        mol = to_neutral(mol)

    if mol:
        # reload molecule, because rd-f***ing-kit
        try:
            # Try catch to avoid occasional aromaticity errors
            # NOTE(hadim): is that still needed?
            return to_mol(to_smiles(mol))
        except Exception:
            return None
    return mol
