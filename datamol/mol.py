from typing import Union
from typing import List
from typing import Tuple
from typing import Optional
from typing import Dict
from typing import Any

import copy
import random

from rdkit import Chem
from rdkit.Chem import rdmolops
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


def copy_mol(mol: Chem.rdchem.Mol) -> Chem.rdchem.Mol:
    """Copy a molecule and return a new one.

    Args:
        mol: a molecule to copy.
    """
    return copy.deepcopy(mol)


def to_mol(
    mol: str,
    add_hs: bool = False,
    explicit_only: bool = False,
    ordered: bool = False,
    kekulize: bool = False,
    sanitize: bool = True,
) -> Optional[Chem.rdchem.Mol]:
    """Convert an input molecule (smiles representation) into a `Chem.rdchem.Mol`.

    Args:
        mol: SMILES of a molecule or a molecule.
        add_hs: Whether hydrogens should be added the molecule.
        explicit_only: Whether to only add explicit hydrogen or both
            (implicit and explicit). when `add_hs` is set to True.
        ordered: Whether the atom should be ordered. This option is
            important if you want to ensure that the features returned will always maintain
            a single atom order for the same molecule, regardless of its original SMILES representation.
        kekulize: Whether to perform kekulization of the input molecules.
        sanitize: Whether to apply rdkit sanitization when input is a SMILES.

    Returns:
        mol: the molecule if some conversion have been made. If the conversion fails
        None is returned so make sure that you handle this case on your own.
    """

    if not isinstance(mol, (str, Chem.rdchem.Mol)):
        raise ValueError(f"Input should be a Chem.rdchem.Mol or a string instead of '{type(mol)}'")

    if isinstance(mol, str):
        _mol = Chem.MolFromSmiles(mol, sanitize=sanitize)

        if not sanitize and _mol is not None:
            _mol.UpdatePropertyCache(False)
    else:
        _mol = mol

    # Add hydrogens
    if _mol is not None and add_hs:
        _mol = Chem.AddHs(_mol, explicitOnly=explicit_only)

    # Reorder atoms
    if _mol is not None and ordered:
        _mol = reorder_atoms(_mol)

    if _mol is not None and kekulize:
        Chem.Kekulize(_mol, clearAromaticFlags=False)
    return _mol


def reorder_atoms(
    mol: Chem.rdchem.Mol,
    break_ties: bool = True,
    include_chirality: bool = True,
    include_isotopes: bool = True,
) -> Optional[Chem.rdchem.Mol]:
    """Reorder the atoms in a mol. It ensures a single atom order for the same molecule,
    regardless of its original representation.

    Args:
        mol: a molecule.
        break_ties: Force breaking of ranked ties.
        include_chirality: Use chiral information when computing rank.
        include_isotopes: Use isotope information when computing rank.

    Returns:
        mol: a molecule.
    """
    if mol.GetNumAtoms() == 0:
        return mol

    new_order = Chem.CanonicalRankAtoms(
        mol,
        breakTies=break_ties,
        includeChirality=include_chirality,
        includeIsotopes=include_isotopes,
    )
    new_order = sorted([(y, x) for x, y in enumerate(new_order)])
    return Chem.RenumberAtoms(mol, [y for (x, y) in new_order])


def randomize_atoms(mol: Chem.rdchem.Mol) -> Optional[Chem.rdchem.Mol]:
    """Randomize the position of the atoms in a mol.

    Args:
        mol: a molecule.

    Returns:
        mol: a molecule.
    """
    if mol.GetNumAtoms() == 0:
        return mol

    atom_indices = list(range(mol.GetNumAtoms()))
    random.shuffle(atom_indices)
    return Chem.RenumberAtoms(mol, atom_indices)


def to_neutral(mol: Chem.rdchem.Mol) -> Optional[Chem.rdchem.Mol]:
    """Neutralize the charge of a molecule.

    Args:
        mol: a molecule.

    Returns:
        mol: a molecule.
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


def sanitize_mol(
    mol: Chem.rdchem.Mol, charge_neutral: bool = False, sanifix: bool = True
) -> Optional[Chem.rdchem.Mol]:
    """Sanitize molecule and fix common errors.

    Warning:
        The procedure includes a SMILES conversion to avoid accasional aromaticity
        errors. In consequence, all the properties and the conformers will be lost.

    Args:
        mol: a molecule.
        charge_neutral: whether charge neutralization should be applied.
        sanifix: whether to run the sanifix from James Davidson
            (sanifix4.py) that try to adjust aromatic nitrogens.

    Returns:
        mol: a molecule.
    """
    if mol is None:
        return mol

    if charge_neutral:
        mol = to_neutral(mol)

    if sanifix:
        mol = _sanifix4.sanifix(mol)

    if mol:
        try:
            # Try catch to avoid occasional aromaticity errors
            return to_mol(dm.to_smiles(mol), sanitize=True)  # type: ignore
        except Exception:
            return None
    return mol


def sanitize_smiles(smiles: str, isomeric: bool = True) -> Optional[str]:
    """Takes SMILES string and returns its sanitized version.

    Args:
        smiles: smiles to be sanitized.
        isomeric: Whether to include information about stereochemistry in the SMILES.

    Returns:
        sanitized smiles.
    """
    try:
        mol = dm.to_mol(smiles, sanitize=False)
        mol = dm.sanitize_mol(mol, False)
    except Exception:
        return None

    if mol is None:
        return None

    try:
        smiles = dm.to_smiles(mol, isomeric=isomeric)  # type: ignore
    except:
        return None
    return smiles


def sanitize_first(mols: List[Chem.rdchem.Mol], charge_neutral: bool = False, sanifix: bool = True):
    """Sanitize a list of molecules and return the first valid molecule seen in the list.

    Args:
        mols: a list of molecules.
        charge_neutral: whether charge neutralization should be applied.
        sanifix: whether to run the sanifix from James Davidson
            (sanifix4.py) that try to adjust aromatic nitrogens.

    Returns:
        mol: a molecule.
    """
    for mol in mols:
        mol = sanitize_mol(mol, charge_neutral=charge_neutral, sanifix=sanifix)
        if mol:
            return mol
    return None


def standardize_smiles(smiles: str, tautomer: bool = False):
    r"""
    Apply smile standardization procedure. This is a convenient function wrapped arrounf RDKit
    smiles standardizer and tautomeric canonicalization.

    Args:
        smiles: Smiles to standardize
        tautomer: Whether to canonicalize tautomers

    Returns:
        standard_smiles: the standardized smiles
    """

    smiles = rdMolStandardize.StandardizeSmiles(smiles)
    if tautomer:
        smiles = canonicalize_tautomer_smiles(smiles)
    return smiles


def standardize_mol(
    mol: Chem.rdchem.Mol,
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
        mol: The molecule to standardize.
        disconnect_metals: Whether to disconnect the metallic atoms from non-metals
        normalize: Whether to apply normalization (correct functional groups and recombine charges).
        reionize: Whether to apply molecule reionization
        uncharge: Whether to remove all charge from molecule
        stereo: Whether to attempt to assign stereochemistry

    Returns:
        mol: The standardized molecule.
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


def fix_valence_charge(mol: Chem.rdchem.Mol, inplace: bool = False) -> Optional[Chem.rdchem.Mol]:
    """Fix valence issues that are due to incorrect charges.

    Args:
        mol: Input molecule with incorrect valence for some atoms
        inplace: Whether to modify in place or make a copy.

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


def incorrect_valence(a: Union[Chem.rdchem.Mol, Chem.rdchem.Atom], update: bool = False) -> bool:
    """Check if an atom connection is not valid or all the atom of a molecule.

    Args:
        a: atom or molecule to check for valence issue.
        update: Update owning molecule property cache first.

    Returns:
        Whether the input atom valence is correct.
    """
    if isinstance(a, Chem.rdchem.Mol):
        a.UpdatePropertyCache(False)
        vm = rdMolStandardize.RDKitValidation()
        return len(vm.validate(a)) > 0

    if update:
        m = a.GetOwningMol()
        m.UpdatePropertyCache(False)
    return (a.GetImplicitValence() == 0) and (
        a.GetExplicitValence() > max(PERIODIC_TABLE.GetValenceList(a.GetSymbol()))
    )


def decrease_bond(bond: Chem.rdchem.Bond) -> Optional[Union[list, Chem.rdchem.Bond]]:
    """Remove one single bond from the input bond. Note that you should
    first kekulize your molecules and remove non-standard bond.

    Args:
        bond: a bond.
    """
    if bond.GetBondType() == TRIPLE_BOND:
        return DOUBLE_BOND
    if bond.GetBondType() == DOUBLE_BOND:
        return SINGLE_BOND
    if bond.GetBondType() == SINGLE_BOND:
        return None
    return bond


def fix_valence(
    mol, inplace: bool = False, allow_ring_break: bool = False
) -> Optional[Chem.rdchem.Mol]:
    """Identify and try to fix valence issues by removing any supplemental bond
    that should not be in the graph.

    Args:
        mol: input molecule with incorrect valence for some atoms
        inplace: Whether to modify in place or make a copy
        allow_ring_break: Whether bond removal involving ring is allowed.

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
        m = Chem.RemoveHs(
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


def adjust_singleton(mol: Chem.rdchem.Mol) -> Optional[Chem.rdchem.Mol]:
    """Remove all atoms that are essentially disconnected singleton nodes in the molecular graph.
    For example, the chlorine atom and methane fragment will be removed in Cl.[N:1]1=CC(O)=CC2CCCCC12.CC.C",
    but not the ethane fragment.

    Args:
        mol: a molecule.
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


def remove_dummies(mol: Chem.rdchem.Mol, dummy: str = "*") -> Optional[Chem.rdchem.Mol]:
    """Remove dummy atoms from molecules."""
    du = dm.to_mol(dummy)
    out = mol
    try:
        out = Chem.ReplaceSubstructs(mol, du, dm.to_mol("[H]"), True)[0]
        out = Chem.RemoveHs(out)
    except Exception as e:
        out = Chem.DeleteSubstructs(mol, du)
    return out


def fix_mol(
    mol: Chem.rdchem.Mol,
    n_iter: int = 1,
    remove_singleton: bool = False,
    largest_only: bool = False,
    inplace: bool = False,
) -> Optional[Chem.rdchem.Mol]:
    """Fix error in molecule using a greedy approach.

    Args:
        mol: input molecule to fix
        n_iter: Number of valence fix iteration to apply
        remove_singleton: Whether `adjust_singleton` should be applied
        largest_only: Whether only the largest fragment should be kept
        inplace: Whether to return a copy of the mol or perform in place operation

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
    mol: Chem.rdchem.Mol,
    atom: str = "C",
    dummy: str = "*",
    replace_all: bool = True,
) -> Optional[Chem.rdchem.Mol]:
    """Remove dummy atoms from molecules.

    Args:
        mol: molecule with dummies
        atom: replacement atom, default is carbon
        dummy: dummy atom representation
        replace_all: Whether to replace all dummies

    Returns:
        mol: Molecule with dummy replaced
    """
    du = Chem.MolFromSmiles(dummy)
    replacement = Chem.MolFromSmiles(atom)
    out = Chem.ReplaceSubstructs(mol, du, replacement, replaceAll=replace_all)[0]
    return out


def keep_largest_fragment(mol: Chem.rdchem.Mol) -> Optional[Chem.rdchem.Mol]:
    """Only keep largest fragment of each molecule."""
    return max(
        rdmolops.GetMolFrags(mol, asMols=True),
        default=mol,
        key=lambda m: m.GetNumAtoms(),
    )


def is_transition_metal(at: Chem.rdchem.Atom) -> bool:
    """Check if atom is a transition metal.

    Args:
        at: an atom.
    """
    n = at.GetAtomicNum()
    return (n >= 22 and n <= 29) or (n >= 40 and n <= 47) or (n >= 72 and n <= 79)


def set_dative_bonds(
    mol: Chem.rdchem.Mol, from_atoms: Tuple[int, int] = (7, 8)
) -> Optional[Chem.rdchem.Mol]:
    """Replaces some single bonds between metals and atoms with atomic numbers in fromAtoms
    with dative bonds. The replacement is only done if the atom has "too many" bonds.

    Arguments:
        mol: molecule with bond to modify
        from_atoms: List of atoms  (symbol or atomic number) to consider for bond replacement.
            By default, only Nitrogen (7) and Oxygen (8) are considered.

    Returns:
        The modified molecule.
    """
    rwmol = Chem.RWMol(mol)
    rwmol.UpdatePropertyCache(strict=False)

    metals = [at for at in rwmol.GetAtoms() if is_transition_metal(at)]
    for metal in metals:
        for nbr in metal.GetNeighbors():
            if (nbr.GetAtomicNum() in from_atoms or nbr.GetSymbol() in from_atoms) and (
                nbr.GetExplicitValence() > PERIODIC_TABLE.GetDefaultValence(nbr.GetAtomicNum())
                and rwmol.GetBondBetweenAtoms(nbr.GetIdx(), metal.GetIdx()).GetBondType()
                == SINGLE_BOND
            ):
                rwmol.RemoveBond(nbr.GetIdx(), metal.GetIdx())
                rwmol.AddBond(nbr.GetIdx(), metal.GetIdx(), DATIVE_BOND)
    return rwmol


def set_mol_props(
    mol: Chem.rdchem.Mol, props: Dict[str, Any], copy: bool = False
) -> Chem.rdchem.Mol:
    """Set properties to a mol from a dict.

    Args:
        mol: the mol where to copy the props.
        props: the props to copy.
        copy: whether to copy the provided mol

    """

    if copy is True:
        mol = dm.copy_mol(mol)

    for k, v in props.items():
        if isinstance(v, bool):
            mol.SetBoolProp(k, v)
        elif isinstance(v, int):
            mol.SetIntProp(k, v)
        elif isinstance(v, float):
            mol.SetDoubleProp(k, v)
        else:
            mol.SetProp(k, str(v))

    return mol


def copy_mol_props(source: Chem.rdchem.Mol, destination: Chem.rdchem.Mol):
    """Copy properties from one source molecule to another destination
    molecule.

    Args:
        source: a molecule to copy from.
        destination: a molecule to copy to.
    """

    props = source.GetPropsAsDict()
    dm.set_mol_props(destination, props)


def enumerate_tautomers(mol: Chem.rdchem.Mol, n_variants: int = 20):
    """Enumerate the possible tautomers of the current molecule.

    Original source: the `openff-toolkit` lib.

    Args:
        mol: The molecule whose state we should enumerate.
        n_variants: The maximum amount of molecules that should be returned.
    """
    # safety first
    mol = copy_mol(mol)

    enumerator = rdMolStandardize.TautomerEnumerator()
    enumerator.SetMaxTautomers(n_variants)
    tautomers = enumerator.Enumerate(mol)
    return list(tautomers)


def enumerate_stereoisomers(
    mol,
    n_variants: int = 20,
    undefined_only: bool = False,
    rationalise: bool = True,
):
    """Enumerate the stereocenters and bonds of the current molecule.

    Original source: the `openff-toolkit` lib.

    Warning: this function can be computationnaly intensive.

    Args:
        mol: The molecule whose state we should enumerate.
        n_variants: The maximum amount of molecules that should be returned.
        undefined_only: If we should enumerate all stereocenters and bonds or only those
            with undefined stereochemistry.
        rationalise: If we should try to build and rationalise the molecule to ensure it
            can exist.
    """
    from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
    from rdkit.Chem.EnumerateStereoisomers import StereoEnumerationOptions

    # safety first
    mol = copy_mol(mol)

    # in case any bonds/centers are missing stereo chem flag it here
    Chem.AssignStereochemistry(mol, force=False, flagPossibleStereoCenters=True, cleanIt=True)  # type: ignore
    Chem.FindPotentialStereoBonds(mol, cleanIt=True)  # type: ignore

    # set up the options
    stereo_opts = StereoEnumerationOptions(
        tryEmbedding=rationalise,
        onlyUnassigned=undefined_only,
        maxIsomers=n_variants,
    )

    try:
        isomers = tuple(EnumerateStereoisomers(mol, options=stereo_opts))
    except:
        # NOTE(hadim): often got "Stereo atoms should be specified before specifying CIS/TRANS bond stereochemistry"
        # for the ligand of reference (coming from the PDB). Not sure how to handle that.
        isomers = []

    variants = []
    for isomer in isomers:
        # isomer has CIS/TRANS tags so convert back to E/Z
        Chem.SetDoubleBondNeighborDirections(isomer)  # type: ignore
        Chem.AssignStereochemistry(isomer, force=True, cleanIt=True)  # type: ignore
        variants.append(isomer)

    return variants


def atom_indices_to_mol(mol: Chem.rdchem.Mol, copy: bool = False):
    """Add the `molAtomMapNumber` property to each atoms.

    Args:
        mol: a molecule
        copy: Whether to copy the molecule.
    """

    if copy is True:
        mol = copy_mol(mol)

    for atom in mol.GetAtoms():
        atom.SetProp("molAtomMapNumber", str(atom.GetIdx()))
    return mol
