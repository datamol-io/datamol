from typing import Union
from typing import List
from typing import Tuple

import copy
import random

from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.MolStandardize import canonicalize_tautomer_smiles

import datamol as dm
from . import _sanifix4

PERIODIC_TABLE = Chem.rdchem.GetPeriodicTable()
TRIPLE_BOND = Chem.rdchem.BondType.TRIPLE
DOUBLE_BOND = Chem.rdchem.BondType.DOUBLE
SINGLE_BOND = Chem.rdchem.BondType.SINGLE
AROMATIC_BOND = Chem.rdchem.BondType.AROMATIC
DATIVE_BOND = Chem.rdchem.BondType.DATIVE


def to_mol(
    mol: str,
    add_hs: bool = False,
    explicit_only: bool = False,
    ordered: bool = False,
    sanitize: bool = True,
):
    """Convert an input molecule (smiles representation) into a `Chem.Mol`.

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


def sanitize_mol(mol: Chem.Mol, charge_neutral: bool = False, sanifix: bool = True):
    """Sanitize molecule and fix common errors.

    Args:
        mol (Chem.Mol): a molecule.
        charge_neutral (bool, optional): whether charge neutralization should be applied.
        sanifix (bool, optional): whether to run the sanifix from James Davidson
            (sanifix4.py) that try to adjust aromatic nitrogens. Default to True.

    Returns:
        mol (Chem.Mol): a molecule.
    """
    if mol is None:
        return mol

    if charge_neutral:
        mol = to_neutral(mol)

    if sanifix:
        mol = _sanifix4.sanifix(mol)

    if mol:
        # reload molecule, because rd-f***ing-kit
        try:
            # Try catch to avoid occasional aromaticity errors
            # NOTE(hadim): is that still needed?
            return to_mol(dm.to_smiles(mol), sanitize=True)
        except Exception:
            return None
    return mol


def sanitize_smiles(smiles: str, isomeric: bool = True):
    """
    Takes list of SMILES strings and returns list of their sanitized versions.

    Args:
        smiles: iterator
            List of smiles to be converted.
        isomeric: bool, optional
            Whether to include information about stereochemistry in the SMILES.
            (Default value=True)
        kwargs: dict
            Any additional parameter to pass to the smile parser

    Returns:
        new_smiles: list
            list of cleaned smiles. Invalid molecules are returned as None
    """
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        mol = dm.sanitize_mol(mol, False)
    except Exception:
        return None

    if mol is None:
        return None

    try:
        smiles = Chem.MolToSmiles(mol, isomericSmiles=isomeric)
    except:
        return None
    return smiles


def sanitize_best(mols: List[Chem.Mol], charge_neutral: bool = False, sanifix: bool = True):
    """Sanitize a list of molecules and return the first valid molecule seen in the list.

    Args:
        mols (List[Chem.Mol]): a list of molecules.
        charge_neutral (bool, optional): whether charge neutralization should be applied.
        sanifix (bool, optional): whether to run the sanifix from James Davidson
            (sanifix4.py) that try to adjust aromatic nitrogens. Default to True.

    Returns:
        mol (Chem.Mol): a molecule.
    """
    for mol in mols:
        mol = sanitize_mol(mol, charge_neutral=charge_neutral, sanifix=sanifix)
        print(mol)
        if mol:
            return mol
    return None


def standardize_smiles(smiles: str, tautomer: bool = False):
    r"""
    Apply smile standardization procedure. This is a convenient function wrapped arrounf RDKit
    smiles standardizer and tautomeric canonicalization.

    Args:
        smiles (str): Smiles to standardize
        tautomer (bool, optional): Whether to canonicalize tautomers

    Returns:
        standard_smiles (str): the standardized smiles
    """

    smiles = rdMolStandardize.StandardizeSmiles(smiles)
    if tautomer:
        smiles = canonicalize_tautomer_smiles(smiles)
    return smiles


def standardize_mol(
    mol: Chem.Mol,
    disconnect_metals: bool = False,
    normalize: bool = True,
    reionize: bool = True,
    uncharge: bool = False,
    stereo: bool = True,
):
    r"""
    This function returns a standardized version the given molecule, with or without disconnect the metals.
    The process is apply in the order of the argument.

    Arguments:
        mol (Chem.Mol): The molecule to standardize.
        disconnect_metals (bool, optional): Whether to disconnect the metallic atoms from non-metals
            Default to False
        normalizer (bool, optional): Whether to apply normalization (correct functional groups and recombine charges).
            Default to True
        reionize (bool, optional): Whether to apply molecule reionization
            Default to True
        uncharge (bool, optional): Whether to remove all charge from molecule
            Default to False
        stereo (bool, optional): Whether to attempt to assign stereochemistry
            Default to True

    Returns:
        mol (Chem.Mol): The standardized molecule.
    """
    mol = copy.copy(mol)

    if disconnect_metals:
        md = rdMolStandardize.MetalDisconnector()
        mol = md.Disconnect(mol)

    if normalize:
        mol = rdMolStandardize.Normalize(mol)

    if reionize:
        reionizer = rdMolStandardize.Reionizer()
        mol = reionizer.reionize(mol)

    if uncharge:
        uncharger = rdMolStandardize.Uncharger()
        mol = uncharger.uncharge(mol)

    if stereo:
        Chem.AssignStereochemistry(mol, force=False, cleanIt=True)

    return mol


def fix_valence_charge(mol: Chem.Mol, inplace: bool = False):
    """Fix valence issues that are due to incorrect charges.

    Args:
        mol (Chem.Mol): Input molecule with incorrect valence for some atoms
        inplace (bool, optional): Whether to modify in place or make a copy.
            Default to False.

    Returns:
        Fixed molecule via charge correction or original molecule if failed.
    """

    vm = rdMolStandardize.RDKitValidation()

    # Don't fix something that is not broken
    if len(vm.validate(mol)) > 0:

        if not inplace:
            mol = copy.copy(mol)

        mol.UpdatePropertyCache(False)
        for a in mol.GetAtoms():
            n_electron = (
                a.GetImplicitValence()
                + a.GetExplicitValence()
                - dm.PERIODIC_TABLE.GetDefaultValence(a.GetSymbol())
            )
            a.SetFormalCharge(n_electron)

    return mol


def incorrect_valence(a: Union[Chem.Mol, Chem.rdchem.Atom], update: bool = False):
    """Check if an atom connection is not valid or all the atom of a molecule.

    Args:
        a (Chem.rdchem.Atom, Chem.Mol): atom or molecule to check for valence issue.
        update (bool, optional): Update owning molecule property cache first.

    Returns:
        Whether the input atom valence is correct.
    """
    if isinstance(a, Chem.Mol):
        a.UpdatePropertyCache(False)
        vm = rdMolStandardize.RDKitValidation()
        return len(vm.validate(a)) > 0

    if update:
        m = a.GetOwningMol()
        m.UpdatePropertyCache(False)
    return (a.GetImplicitValence() == 0) and (
        a.GetExplicitValence() > max(PERIODIC_TABLE.GetValenceList(a.GetSymbol()))
    )


def decrease_bond(bond: Chem.rdchem.Bond):
    """Remove one single bond from the input bond. Note that you should
    first kekulize your molecules and remove non-standard bond.
    """
    if bond.GetBondType() == TRIPLE_BOND:
        return DOUBLE_BOND
    if bond.GetBondType() == DOUBLE_BOND:
        return SINGLE_BOND
    if bond.GetBondType() == SINGLE_BOND:
        return None
    return bond


def fix_valence(mol, inplace: bool = False, allow_ring_break: bool = False):
    """Identify and try to fix valence issues by removing any supplemental bond
    that should not be in the graph.

    Args:
        mol (Chem.Mol): input molecule with incorrect valence for some atoms
        inplace (bool, optional): Whether to modify in place or make a copy
            Default to False.
        allow_ring_break (bool, optional): Whether bond removal involving ring is allowed.
            Default to False.

    Returns:
        Fixed potential valence issue in molecule or original molecule when nothing is broken
        of if failed.
    """
    if not inplace:
        mol = copy.copy(mol)

    vm = rdMolStandardize.RDKitValidation()
    if len(vm.validate(mol)) == 0:  # don't fix something that is not broken
        return mol

    try:
        m = Chem.rdmolops.RemoveHs(
            mol,
            implicitOnly=False,
            updateExplicitCount=True,
            sanitize=False,
        )
        m.UpdatePropertyCache(False)

        # first pass using explicit false count
        for atom in m.GetAtoms():
            while incorrect_valence(atom) and atom.GetTotalNumHs() > 0:
                cur_hydrogen = atom.GetTotalNumHs()
                atom.SetNumExplicitHs(max(0, cur_hydrogen - 1))
                atom.SetFormalCharge(max(0, atom.GetFormalCharge() - 1))
                # atom.SetNumRadicalElectrons(0)
            atom.UpdatePropertyCache(False)

        em = Chem.RWMol(m)
        bonds = em.GetBonds()
        bonds = [
            bond
            for bond in bonds
            if any(
                [
                    incorrect_valence(bond.GetBeginAtom()),
                    incorrect_valence(bond.GetEndAtom()),
                ]
            )
        ]
        for bond in bonds:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if incorrect_valence(a1) or incorrect_valence(a2):
                mbond = decrease_bond(bond)
                if allow_ring_break or (mbond or not bond.IsInRing()):
                    em.RemoveBond(a1.GetIdx(), a2.GetIdx())
                    if mbond is not None:
                        em.AddBond(a1.GetIdx(), a2.GetIdx(), mbond)
            a1.UpdatePropertyCache(False)
            a2.UpdatePropertyCache(False)
        m = em.GetMol()

    except Exception:
        return None

    return m


def adjust_singleton(mol: Chem.Mol):
    """Remove all atoms that are essentially disconnected singleton nodes in the molecular graph.
    For example, the chlorine atom and methane fragment will be removed in Cl.[N:1]1=CC(O)=CC2CCCCC12.CC.C",
    but not the ethane fragment.
    """
    to_rem = []
    em = Chem.RWMol(mol)
    for atom in mol.GetAtoms():
        if atom.GetExplicitValence() == 0:
            to_rem.append(atom.GetIdx())
    to_rem.sort(reverse=True)
    for a_idx in to_rem:
        em.RemoveAtom(a_idx)
    return em.GetMol()


def remove_dummies(mol: Chem.Mol, dummy: str = "*"):
    """Remove dummy atoms from molecules."""
    du = Chem.MolFromSmiles(dummy)
    out = mol
    try:
        out = Chem.ReplaceSubstructs(mol, du, Chem.MolFromSmiles("[H]"), True)[0]
        out = Chem.RemoveHs(out)
    except Exception as e:
        out = Chem.DeleteSubstructs(mol, du)
    return out


def fix_mol(
    mol: Chem.Mol,
    n_iter: int = 1,
    remove_singleton: bool = False,
    largest_only: bool = False,
    inplace: bool = False,
):
    """Fix error in molecule using a greedy approach.

    Args:
        mol (Chem.Mol): input molecule to fix
        n_iter (int, optional): Number of valence fix iteration to apply
            Default to 1.
        remove_singleton (bool, optional): Whether `adjust_singleton` should be applied
            Default to False.
        largest_only (bool, optional): Whether only the largest fragment should be kept
            Default to False.
        inplace (bool, optional): Whether to return a copy of the mol or perform in place operation
            Default to False.

    Returns:
        Fixed molecule.
    """

    if not inplace:
        mol = copy.copy(mol)

    m = sanitize_mol(mol) or mol  # fail back to mol when the fixer fail

    if m is not None:
        m = remove_dummies(m)
        for _ in range(n_iter):
            m = fix_valence(m)

        if remove_singleton:
            m = adjust_singleton(m)

        if largest_only:
            # m = max(Chem.rdmolops.GetMolFrags(m, asMols=True, sanitizeFrags=False), key=lambda m: m.GetNumAtoms())
            m = rdMolStandardize.FragmentParent(m, skipStandardize=True)

    return m


def replace_dummies_atoms(
    mol: Chem.Mol,
    atom: str = "C",
    dummy: str = "*",
    replace_all: bool = True,
):
    r"""Remove dummy atoms from molecules.

    Args:
        mol: <Chem.Mol>): molecule with dummies
        atom (str, optional): replacement atom, default is carbon
            Default to'C'
        dummy (str, optional): dummy atom representation
            Default to '*'.
        replace_all (bool, optional): Whether to replace all dummies

    Returns:
        mol (Chem.Mol): Molecule with dummy replaced
    """
    du = Chem.MolFromSmiles(dummy)
    replacement = Chem.MolFromSmiles(atom)
    out = Chem.ReplaceSubstructs(mol, du, replacement, replaceAll=replace_all)[0]
    return out


def keep_largest_fragment(mol: Chem.Mol):
    """Only keep largest fragment of each molecule."""
    return max(
        rdmolops.GetMolFrags(mol, asMols=True),
        default=mol,
        key=lambda m: m.GetNumAtoms(),
    )


def is_transition_metal(at: Chem.rdchem.Atom):
    """Check if atom is a transition metal."""
    n = at.GetAtomicNum()
    return (n >= 22 and n <= 29) or (n >= 40 and n <= 47) or (n >= 72 and n <= 79)


def set_dative_bonds(mol: Chem.Mol, from_atoms: Tuple[int, int] = (7, 8)):
    """Replaces some single bonds between metals and atoms with atomic numbers in fromAtoms
    with dative bonds. The replacement is only done if the atom has "too many" bonds.

    Arguments:
        mol (Chem.Mol): molecule with bond to modify
        from_atoms (list or tuple): List of atoms  (symbol or atomic number) to consider for bond replacement.
            By default, only Nitrogen (7) and Oxygen (8) are considered.

    Returns:
        The modified molecule.
    """
    pt = Chem.GetPeriodicTable()
    rwmol = Chem.RWMol(mol)
    rwmol.UpdatePropertyCache(strict=False)

    metals = [at for at in rwmol.GetAtoms() if is_transition_metal(at)]
    for metal in metals:
        for nbr in metal.GetNeighbors():
            if (nbr.GetAtomicNum() in from_atoms or nbr.GetSymbol() in from_atoms) and (
                nbr.GetExplicitValence() > pt.GetDefaultValence(nbr.GetAtomicNum())
                and rwmol.GetBondBetweenAtoms(nbr.GetIdx(), metal.GetIdx()).GetBondType()
                == SINGLE_BOND
            ):
                rwmol.RemoveBond(nbr.GetIdx(), metal.GetIdx())
                rwmol.AddBond(nbr.GetIdx(), metal.GetIdx(), DATIVE_BOND)
    return rwmol
