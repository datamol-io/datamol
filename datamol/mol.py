from typing import Union

import random
import re

from rdkit import Chem
import selfies as sf


PERIODIC_TABLE = Chem.rdchem.GetPeriodicTable()


def to_mol(
    mol: str,
    add_hs: bool = False,
    explicit_only: bool = False,
    ordered: bool = False,
    sanitize: bool = True,
):
    """Convert an input molecule (smiles representation) into a `Chem.Mol`.

    NOTE(hadim): should we support SELFIES here (deepsmiles)?

    Args:
        mol (str): SMILES of a molecule or a molecule.
        add_hs (bool, optional): Whether hydrogens should be added the molecule. Default to False.
        explicit_only (bool, optional): Whether to only add explicit hydrogen or both
            (implicit and explicit). when `add_hs` is set to True. Default to False.
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


def randomize_atoms(mol: Chem.Mol):
    """Randomize the position of the atoms in a mol.

    Args:
        mol (Chem.Mol): a molecule.

    Returns:
        mol (Chem.Mol): a molecule.
    """
    if mol.GetNumAtoms() == 0:
        return mol

    atom_indices = list(range(mol.GetNumAtoms()))
    random.shuffle(atom_indices)
    return Chem.RenumberAtoms(mol, atom_indices)


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


def sanitize_mol(mol: Chem.Mol, charge_neutral: bool = False):
    """Sanitize molecule and fix common errors.

    Args:
        mol (Chem.Mol): a molecule.
        charge_neutral (bool, optional): whether charge neutralization should be applied.

    Returns:
        mol (Chem.Mol): a molecule.
    """
    if charge_neutral:
        mol = to_neutral(mol)

    if mol:
        # reload molecule, because rd-f***ing-kit
        try:
            # Try catch to avoid occasional aromaticity errors
            # NOTE(hadim): is that still needed?
            return to_mol(to_smiles(mol), sanitize=True)
        except Exception:
            return None
    return mol


def to_smiles(
    mol: Chem.Mol,
    canonical: bool = True,
    isomeric: bool = True,
    ordered: bool = False,
    explicit_bonds: bool = False,
    explicit_hs: bool = False,
    randomize: bool = False,
):
    """Convert a mol to a SMILES.

    Args:
        mol (Chem.Mol): a molecule.
        add_hs (bool, optional): Whether hydrogens should be added the SMILES. Default to False.
        canonical (bool, optional): if false no attempt will be made to canonicalize the molecule.
            Defaults to true.
        isomeric (bool, optional): whether to include information about stereochemistry in
            the SMILES. Defaults to True.
        ordered (bool, optional): whether to force reordering of the atoms
            first. Defaults to False.
        explicit_bonds (bool, optional): if true, all bond orders will be explicitly indicated in
            the output SMILES. Defaults to false.
        explicit_hs (bool, optional): if true, all H counts will be explicitly indicated in the
            output SMILES. Defaults to false.
        randomize (bool, optional): whether to randomize the generated smiles. Override `canonical`.
            Defaults to false.
    """
    if ordered:
        mol = reorder_atoms(mol)

    if randomize:
        mol = randomize_atoms(mol)
        canonical = False

    smiles = None
    try:
        smiles = Chem.MolToSmiles(
            mol,
            isomericSmiles=isomeric,
            canonical=canonical,
            allBondsExplicit=explicit_bonds,
            allHsExplicit=explicit_hs,
        )
    except:
        return None
    return smiles


def to_selfies(mol: Union[str, Chem.Mol]):
    """Convert a mol to SELFIES.

    Args:
        mol (Chem.Mol or str): a molecule or a SMILES.

    Returns:
        selfies (str): SELFIES string.
    """
    if mol is None:
        return None

    if isinstance(mol, Chem.Mol):
        mol = to_smiles(mol)

    return sf.encoder(mol)


def from_selfies(selfies: str, as_mol: bool = False):
    """Convert a SEFLIES to a smiles or a mol.

    Args:
        selfies (str): a selfies.
        as_mol (str, optional): whether to return a mol or a smiles.

    Returns:
        smiles or mol (str, Chem.Mol))
    """
    if selfies is None:
        return None

    smiles = sf.decoder(selfies)

    if as_mol and smiles is not None:
        return to_mol(smiles)

    return smiles


def to_smarts(mol: Union[str, Chem.Mol], keep_hs: bool = True):
    """Convert a molecule to a smarts.

    Args:
        mol (Chem.Mol): a molecule.
        keep_hs (bool, optional): Whether to keep hydrogen. This will increase the count of H atoms
            for atoms with attached hydrogens to create a valid smarts.
            e.g. [H]-[CH2]-[*] -> [H]-[CH3]-[*]

    Returns:
        smarts of the molecule
    """

    if mol is None:
        return None

    # Change the isotope to 42
    for atom in mol.GetAtoms():
        if keep_hs:
            s = sum(na.GetAtomicNum() == 1 for na in atom.GetNeighbors())
            if s:
                atom.SetNumExplicitHs(atom.GetTotalNumHs() + s)
        atom.SetIsotope(42)

    # Print out the smiles, all the atom attributes will be fully specified
    smarts = to_smiles(mol, isomeric=True, explicit_bonds=True)

    if smarts is None:
        return None

    # Remove the 42 isotope labels
    smarts = re.sub(r"\[42", "[", smarts)
    return smarts
