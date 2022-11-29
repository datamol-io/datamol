from typing import Union
from typing import List
from typing import Tuple
from typing import Optional
from typing import Dict
from typing import Any
from typing import Set

import copy
import random
import itertools
import hashlib

from loguru import logger

import numpy as np
import numpy.typing as npt
from rdkit import Chem

from rdkit.Chem import rdmolops
from rdkit.Chem import rdchem
from rdkit.Chem import rdmolfiles
from rdkit.Chem import rdGeometry

from rdkit.Chem.Scaffolds import MurckoScaffold

from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.MolStandardize import canonicalize_tautomer_smiles

from . import _sanifix4
from .types import Mol
from .convert import to_inchikey_non_standard
from .convert import to_inchikey
from .convert import to_smiles
from .convert import from_smarts
from .log import without_rdkit_log
from ._version import is_lower_than_current_rdkit_version


PERIODIC_TABLE = Chem.rdchem.GetPeriodicTable()
TRIPLE_BOND = Chem.rdchem.BondType.TRIPLE
DOUBLE_BOND = Chem.rdchem.BondType.DOUBLE
SINGLE_BOND = Chem.rdchem.BondType.SINGLE
AROMATIC_BOND = Chem.rdchem.BondType.AROMATIC
DATIVE_BOND = Chem.rdchem.BondType.DATIVE
UNSPECIFIED_BOND = Chem.rdchem.BondType.UNSPECIFIED


def copy_mol(mol: Mol) -> Mol:
    """Copy a molecule and return a new one.

    Args:
        mol: a molecule to copy.
    """
    return copy.deepcopy(mol)


def to_mol(
    mol: Union[str, Mol],
    add_hs: bool = False,
    explicit_only: bool = False,
    ordered: bool = False,
    kekulize: bool = False,
    sanitize: bool = True,
    allow_cxsmiles: bool = True,
    parse_name: bool = True,
    remove_hs: bool = True,
    strict_cxsmiles: bool = True,
) -> Optional[Mol]:
    """Convert an input molecule (smiles representation) into a `Mol`.

    Args:
        mol: A SMILES or a molecule.
        add_hs: Whether hydrogens should be added the molecule after the SMILES has been parsed.
        explicit_only: Whether to only add explicit hydrogen or both
            (implicit and explicit). when `add_hs` is set to True.
        ordered: Whether the atom should be ordered. This option is
            important if you want to ensure that the features returned will always maintain
            a single atom order for the same molecule, regardless of its original SMILES representation.
        kekulize: Whether to perform kekulization of the input molecules.
        sanitize: Whether to apply rdkit sanitization when input is a SMILES.
        allow_cxsmiles: Recognize and parse CXSMILES.
        parse_name: Parse (and set) the molecule name as well.
        remove_hs: Wether to remove the hydrogens from the input SMILES.
        strict_cxsmiles: Throw an exception if the CXSMILES parsing fails.

    Returns:
        mol: the molecule if some conversion have been made. If the conversion fails
        None is returned so make sure that you handle this case on your own.
    """

    if not isinstance(mol, (str, Mol)):
        raise ValueError(f"Input should be a Mol or a string instead of '{type(mol)}'")

    if isinstance(mol, str):
        smiles_params = rdmolfiles.SmilesParserParams()
        smiles_params.sanitize = sanitize
        smiles_params.allowCXSMILES = allow_cxsmiles
        smiles_params.parseName = parse_name
        smiles_params.removeHs = remove_hs
        smiles_params.strictCXSMILES = strict_cxsmiles

        _mol = rdmolfiles.MolFromSmiles(mol, params=smiles_params)

        if not sanitize and _mol is not None:
            _mol.UpdatePropertyCache(False)
    else:
        _mol = mol

    # Add hydrogens
    if _mol is not None and add_hs:
        _mol = rdmolops.AddHs(_mol, explicitOnly=explicit_only, addCoords=True)

    # Reorder atoms
    if _mol is not None and ordered:
        _mol = reorder_atoms(_mol)

    if _mol is not None and kekulize:
        rdmolops.Kekulize(_mol, clearAromaticFlags=False)

    return _mol


def same_mol(
    mol1: Optional[Mol],
    mol2: Optional[Mol],
    use_non_standard_inchikey: bool = False,
) -> bool:
    """Check two molecules are the same by comparing their InChiKey.

    Invalid molecules (None) are always considered as not the same.

    Args:
        mol1: A molecule.
        mol2: A molecule.
        use_non_standard_inchikey: Whether to use the standard or non-standard InChiKey.
    """

    if mol1 is None or mol2 is None:
        return False

    if use_non_standard_inchikey:
        return to_inchikey_non_standard(mol1) == to_inchikey_non_standard(mol2)
    else:
        return to_inchikey(mol1) == to_inchikey(mol2)


def unique_id(mol: Mol) -> Optional[str]:
    """A datamol unique molecule ID.

    The ID is an MD5 hash of the non-standard InChiKey provided
    by `dm.to_inchikey_non_standard()`. It guarantees uniqueness for
    different tautomeric forms of the same molecule.

    Args:
        mol: A molecule.
    """
    ik = to_inchikey_non_standard(mol)

    if ik is None:
        return None

    return hashlib.md5(ik.encode("utf-8")).hexdigest()


def hash_mol(mol: Mol, hash_scheme: str = "all") -> str:
    """Generate a unique hash code for a molecule based on chemistry. If two molecules are
    chemically “the same”, they should have the same hash.

    Using molhash adds value beyond using SMILES because it:

    - Ignores SMILES features that are not chemically meaningful (e.g. atom map numbers).
    - Canonicalizes enhanced stereochemistry groups. For example `C[C@H](O)CC |&1:1|` and `C[C@@H](O)CC |&1:1|`
    have the same molhash.
    - Canonicalizes S group data (for example, polymer data).

    Args:
        mol: A molecule.
        hash_scheme: There are three hash schemes:
            - "all": most strict hash scheme utilizing all layers.
            - "no_stereo": excludes stereo sensitive layers.
            - "no_tautomers": excludes tautomer sensitive layers.

    Returns:
        The hash as a 40 chars string.
    """

    if is_lower_than_current_rdkit_version("2022.09"):
        raise NotImplementedError("`datamol.hash_mol() is only available with RDKit>=2022.09.")

    from rdkit.Chem import RegistrationHash

    # Compute the mol layers
    all_layers = RegistrationHash.GetMolLayers(mol)

    # Select the hash scheme
    if hash_scheme == "all":
        hash_scheme_enum = RegistrationHash.HashScheme.ALL_LAYERS
    elif hash_scheme == "no_stereo":
        hash_scheme_enum = RegistrationHash.HashScheme.STEREO_INSENSITIVE_LAYERS
    elif hash_scheme == "no_tautomers":
        hash_scheme_enum = RegistrationHash.HashScheme.TAUTOMER_INSENSITIVE_LAYERS
    else:
        raise ValueError(
            f"`hash_scheme` is invalid. Please choose from: 'all', 'no_stereo' or 'no_tautomers'."
        )

    # Generate the hash
    return RegistrationHash.GetMolHash(all_layers=all_layers, hash_scheme=hash_scheme_enum)


def reorder_atoms(
    mol: Mol,
    break_ties: bool = True,
    include_chirality: bool = True,
    include_isotopes: bool = True,
) -> Optional[Mol]:
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

    new_order = rdmolfiles.CanonicalRankAtoms(
        mol,
        breakTies=break_ties,
        includeChirality=include_chirality,
        includeIsotopes=include_isotopes,
    )
    new_order = sorted([(y, x) for x, y in enumerate(new_order)])
    return rdmolops.RenumberAtoms(mol, [y for (x, y) in new_order])


def randomize_atoms(mol: Mol) -> Optional[Mol]:
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
    return rdmolops.RenumberAtoms(mol, atom_indices)


def to_neutral(mol: Optional[Mol]) -> Optional[Mol]:
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
    mol: Mol,
    charge_neutral: bool = False,
    sanifix: bool = True,
    verbose: bool = True,
    add_hs: bool = False,
) -> Optional[Mol]:
    """An augmented version of RDKit `sanitize=True`. It uses a
    mol-SMILES-mol conversion to catch potential aromaticity errors
    and try to fix aromatic nitrogen (using the popular sanifix4 script).
    Optionally, it can neutralize the charge of the molecule.

    Note #1: Only the first conformer (if present) will be preserved and
    a warning will be displayed if more than one conformer is detected.

    Note #2: The molecule's properties will be preserved but the atom's
    properties will be lost.

    Args:
        mol: a molecule.
        charge_neutral: whether charge neutralization should be applied.
        sanifix: whether to run the sanifix from James Davidson
            (sanifix4.py) that try to adjust aromatic nitrogens.
        verbose: Whether displaying a warning about multiple conformers.
        add_hs: Add hydrogens to the returned molecule. Useful when the input
            molecule already contains hydrogens.

    Returns:
        mol: a molecule.
    """
    if mol is None:
        return mol

    # Extract properties.
    original_mol = copy_mol(mol)
    properties = original_mol.GetPropsAsDict()

    if charge_neutral:
        mol = to_neutral(mol)

    if sanifix:
        mol = _sanifix4.sanifix(mol)

    if mol is not None:

        # Detect multiple conformers
        if verbose and mol.GetNumConformers() > 1:
            logger.warning(
                f"The molecule contains multiple conformers. Only the first one will be preserved."
            )

        # Try catch to avoid occasional aromaticity errors
        try:
            # `cxsmiles` is used here to preserve the first conformer.
            mol = to_mol(to_smiles(mol, cxsmiles=True), sanitize=True, add_hs=add_hs)
        except Exception:
            mol = None

    if mol is not None:
        # Insert back properties.
        mol = set_mol_props(mol, properties)

    return mol


def sanitize_smiles(smiles: Optional[str], isomeric: bool = True) -> Optional[str]:
    """Takes SMILES string and returns its sanitized version.

    Args:
        smiles: smiles to be sanitized.
        isomeric: Whether to include information about stereochemistry in the SMILES.

    Returns:
        sanitized smiles.
    """

    mol = None

    try:
        mol = to_mol(smiles, sanitize=False)
        mol = sanitize_mol(mol, False)
    except Exception:
        return None

    if mol is None:
        return None

    try:
        smiles = to_smiles(mol, isomeric=isomeric)
    except:
        return None

    return smiles


def sanitize_first(mols: List[Mol], charge_neutral: bool = False, sanifix: bool = True) -> Mol:
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


def standardize_smiles(smiles: str, tautomer: bool = False) -> str:
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
    mol: Mol,
    disconnect_metals: bool = False,
    normalize: bool = True,
    reionize: bool = True,
    uncharge: bool = False,
    stereo: bool = True,
) -> Mol:
    r"""
    This function returns a standardized version the given molecule. It relies on the
    RDKit [`rdMolStandardize` module](https://www.rdkit.org/docs/source/rdkit.Chem.MolStandardize.rdMolStandardize.html)
    which is largely inspired from [MolVS](https://github.com/mcs07/MolVS).

    Arguments:
        mol: A molecule to standardize.

        disconnect_metals: Disconnect metals that are defined as covalently bonded to non-metal.
            Depending on the source of the database, some compounds may be reported in salt form
            or associated to metallic ions (e.g. the sodium salt of a carboxylic compound).
            In most cases, these counter-ions are not relevant so the use of this function is required
            before further utilization of the dataset. In summary the process is the following:

            - Break covalent bonds between metals and organic atoms under certain conditions.
            - First, disconnect N, O, F from any metal. Then disconnect other non-metals from transition metals (with exceptions).
            - For every bond broken, adjust the charges of the begin and end atoms accordingly.

        normalize: Applies a series of standard transformations to correct functional groups and recombine charges.
            It corrects drawing errors and standardizes functional groups in the molecule as well as ensuring the
            overall proper charge of the compound. It includes:

            - Uncharge-separate sulfones
            - Charge-separate nitro groups
            - Charge-separate pyridine oxide
            - Charge-separate azide
            - Charge-separate diazo and azo groups
            - Charge-separate sulfoxides
            - Hydrazine-diazonium system

        reionize: If one or more acidic functionalities are present in the molecule, this option ensures the correct
            neutral/ionized state for such functional groups. Molecules are uncharged by adding and/or removing hydrogens.
            For zwitterions, hydrogens are moved to eliminate charges where possible. However, in cases where there is a
            positive charge that is not neutralizable, an attempt is made to also preserve the corresponding negative charge
            The algorithm works as follows:

            - Use SMARTS to find the strongest protonated acid and the weakest ionized acid.
            - If the ionized acid is weaker than the protonated acid, swap proton and repeat.

        uncharge: This option neutralize the molecule by reversing the protonation state of protonated and deprotonated groups,
            if present (e.g. a carboxylate is re-protonated to the corresponding carboxylic acid).
            In cases where there is a positive charge that is not neutralizable, an attempt is made to also preserve the
            corresponding negative charge to ensure a net zero charge.

        stereo: Stereochemical information is corrected and/or added if missing using built-in RDKit functionality to force a clean recalculation of stereochemistry (`AssignStereochemistry`).

    Returns:
        mol: A standardized molecule.
    """
    mol = copy_mol(mol)

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
        rdmolops.AssignStereochemistry(mol, force=False, cleanIt=True)

    return mol


def fix_valence_charge(mol: Mol, inplace: bool = False) -> Optional[Mol]:
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
                - PERIODIC_TABLE.GetDefaultValence(a.GetSymbol())
            )
            a.SetFormalCharge(n_electron)

    return mol


def incorrect_valence(a: Union[Mol, Chem.rdchem.Atom], update: bool = False) -> bool:
    """Check if an atom connection is not valid or all the atom of a molecule.

    Args:
        a: atom or molecule to check for valence issue.
        update: Update owning molecule property cache first.

    Returns:
        Whether the input atom valence is correct.
    """
    if isinstance(a, Mol):
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


def fix_valence(mol: Mol, inplace: bool = False, allow_ring_break: bool = False) -> Optional[Mol]:
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
        m = remove_hs(
            mol,
            implicit_only=False,
            update_explicit_count=True,
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

        em = rdchem.RWMol(m)
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


def adjust_singleton(mol: Mol) -> Optional[Mol]:
    """Remove all atoms that are essentially disconnected singleton nodes in the molecular graph.
    For example, the chlorine atom and methane fragment will be removed in Cl.[N:1]1=CC(O)=CC2CCCCC12.CC.C",
    but not the ethane fragment.

    Args:
        mol: a molecule.
    """
    to_rem = []
    em = rdchem.RWMol(mol)
    for atom in mol.GetAtoms():
        if atom.GetExplicitValence() == 0:
            to_rem.append(atom.GetIdx())
    to_rem.sort(reverse=True)
    for a_idx in to_rem:
        em.RemoveAtom(a_idx)
    return em.GetMol()


def remove_dummies(mol: Mol, dummy: str = "*") -> Optional[Mol]:
    """Remove dummy atoms from molecules."""

    du = to_mol(dummy)
    out = mol

    try:
        out = rdmolops.ReplaceSubstructs(mol, du, to_mol("[H]"), True)[0]
        out = remove_hs(out)
    except Exception:
        out = rdmolops.DeleteSubstructs(mol, du)
    return out


def fix_mol(
    mol: Mol,
    n_iter: int = 1,
    remove_singleton: bool = False,
    largest_only: bool = False,
    inplace: bool = False,
) -> Optional[Mol]:
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
    mol: Mol,
    atom: str = "C",
    dummy: str = "*",
    replace_all: bool = True,
) -> Optional[Mol]:
    """Remove dummy atoms from molecules.

    Args:
        mol: molecule with dummies
        atom: replacement atom, default is carbon
        dummy: dummy atom representation
        replace_all: Whether to replace all dummies

    Returns:
        mol: Molecule with dummy replaced
    """
    du = to_mol(dummy)
    replacement = to_mol(atom)
    out = rdmolops.ReplaceSubstructs(mol, du, replacement, replaceAll=replace_all)[0]
    return out


def keep_largest_fragment(mol: Mol) -> Optional[Mol]:
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


def set_dative_bonds(mol: Mol, from_atoms: Tuple[int, int] = (7, 8)) -> Optional[Mol]:
    """Replaces some single bonds between metals and atoms with atomic numbers in fromAtoms
    with dative bonds. The replacement is only done if the atom has "too many" bonds.

    Arguments:
        mol: molecule with bond to modify
        from_atoms: List of atoms  (symbol or atomic number) to consider for bond replacement.
            By default, only Nitrogen (7) and Oxygen (8) are considered.

    Returns:
        The modified molecule.
    """
    rwmol = rdchem.RWMol(mol)
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
    mol: Mol,
    props: Dict[str, Any],
    copy: bool = False,
) -> Mol:
    """Set properties to a mol from a dict.

    Args:
        mol: the mol where to copy the props.
        props: the props to copy.
        copy: whether to copy the provided mol

    """

    if copy is True:
        mol = copy_mol(mol)

    for k, v in props.items():
        if isinstance(v, bool):
            mol.SetBoolProp(k, v)
        elif isinstance(v, int):
            # NOTE(hadim): A Python integer is 32 bits and RDKit seems
            # to overflow before that. Here we catch the error
            # and instead uses silently `SetDoubleProp` instead.
            try:
                mol.SetIntProp(k, v)
            except OverflowError:
                mol.SetDoubleProp(k, v)
        elif isinstance(v, float):
            mol.SetDoubleProp(k, v)
        else:
            mol.SetProp(k, str(v))

    return mol


def copy_mol_props(
    source: Mol,
    destination: Mol,
    include_private: bool = False,
    include_computed: bool = False,
):
    """Copy properties from one source molecule to another destination
    molecule.

    Args:
        source: a molecule to copy from.
        destination: a molecule to copy to.
        include_private: Include private properties.
        include_computed: Include computed properties.
    """

    props = source.GetPropsAsDict(includePrivate=include_private, includeComputed=include_computed)
    set_mol_props(destination, props)


def clear_mol_props(
    mol: Mol,
    property_keys: Optional[Union[List[str], str]] = None,
    copy: bool = True,
    include_private: bool = False,
    include_computed: bool = False,
) -> Mol:
    """Clear properties from a molecule.

    Args:
        mol: A molecule.
        property_keys: If set, only the set properties will be cleared.
            It not set, all the properties are cleared.
        copy: Whether to copy the molecule.
        include_private: Whether to also clean the private properties.
        include_computed: Whether to also clean the computed properties.
    """

    if copy:
        mol = copy_mol(mol)

    if property_keys is not None and isinstance(property_keys, str):
        property_keys = [property_keys]

    if property_keys is None:
        props = mol.GetPropsAsDict(includePrivate=include_private, includeComputed=include_computed)
        property_keys = list(props.keys())

    for key in property_keys:
        mol.ClearProp(key)

    return mol


def clear_atom_props(
    mol: Mol,
    property_keys: Optional[Union[List[str], str]] = None,
    copy: bool = True,
    include_private: bool = False,
    include_computed: bool = False,
) -> Mol:
    """Clear atom properties from a molecule.

    Args:
        mol: A molecule.
        property_keys: If set, only the set properties will be cleared.
            It not set, all the properties are cleared.
        copy: Whether to copy the molecule.
        include_private: Whether to also clean the private properties.
        include_computed: Whether to also clean the computed properties.
    """

    if copy:
        mol = copy_mol(mol)

    if property_keys is not None and isinstance(property_keys, str):
        property_keys = [property_keys]

    for atom in mol.GetAtoms():

        if property_keys is None:
            props = atom.GetPropsAsDict(
                includePrivate=include_private, includeComputed=include_computed
            )
            property_keys = list(props.keys())

        for key in property_keys:
            atom.ClearProp(key)

    return mol


def clear_atom_map_number(mol: Mol, copy: bool = True):
    """Clear the `molAtomMapNumber` property of the atom's molecule.

    Args:
        mol: A molecule.
        copy: Whether to copy the molecule.
    """
    return clear_atom_props(mol, copy=copy, property_keys="molAtomMapNumber")


def atom_indices_to_mol(mol: Mol, copy: bool = True):
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


def atom_list_to_bond(
    mol: Mol,
    atom_indices: List[int],
    bond_as_idx: bool = False,
):
    """Return a list of existing bond indices between a list of
    atom indices.

    Args:
        mol: A molecule.
        atom_indices: A list of atom indices.
    """

    # Build an atom map
    atom_map = {}
    submol = rdmolops.PathToSubmol(mol, atom_indices, useQuery=True, atomMap=atom_map)
    atom_map_reversed = {v: k for k, v in atom_map.items()}

    bonds = []

    for bond in submol.GetBonds():
        a1, a2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        ori_a1 = atom_map_reversed[a1]
        ori_a2 = atom_map_reversed[a2]

        if ori_a1 in atom_indices and ori_a2 in atom_indices:
            ori_bond = mol.GetBondBetweenAtoms(ori_a1, ori_a2)
            if bond_as_idx:
                bonds.append(ori_bond.GetIdx())
            else:
                bonds.append(ori_bond)

    return bonds


def substructure_matching_bonds(mol: Mol, query: Mol, **kwargs: Any) -> Tuple[list, list]:
    """Perform a substructure match using `GetSubstructMatches` but instead
    of returning only the atom indices also return the bond indices.

    Args:
        mol: A molecule.
        query: A molecule used as a query to match against.
        **kwargs: Any other arguments to pass to `mol.GetSubstructMatches()`.

    Returns:
        atom_matches: A list of lists of atom indices.
        bond_matches: A list of lists of bond indices.
    """

    # NOTE(hadim): If more substructure functions are added here, consider moving it to
    # a dedicated `substructure` module.

    # Set default arguments
    kwargs.setdefault("uniquify", True)

    # Get the matching atom indices
    atom_matches = list(mol.GetSubstructMatches(query, **kwargs))

    # Get the bond to highligh from the query
    query_bond_indices = [
        (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) for bond in query.GetBonds()
    ]

    # Retrieve the atom indices
    query_atom_indices = [atom.GetIdx() for i, atom in enumerate(query.GetAtoms())]

    bond_matches = []

    for match in atom_matches:

        # Map the atom of the query to the atom of the mol matching the query
        atom_map = dict(zip(query_atom_indices, match))

        # For this match atoms we now, we use the map to retrieve the matching bonds
        # in the mol.
        mol_bond_indices = [(atom_map[a1], atom_map[a2]) for a1, a2 in query_bond_indices]

        # Convert the bond atom indices to bond indices
        mol_bond_indices = [mol.GetBondBetweenAtoms(a1, a2).GetIdx() for a1, a2 in mol_bond_indices]

        bond_matches.append(mol_bond_indices)

    return atom_matches, bond_matches


def protect_atoms(
    mol: Mol,
    substruct: Optional[Mol] = None,
    atoms: Optional[Union[List[int], int]] = None,
    in_place: bool = False,
) -> Mol:
    """Protect a list of atoms or substruct in a molecule.

    The _protected attributes of a molecule is used by RDKit in several functions, especially for reactions
    where "protected" atoms are disallowed from taking part in reactions.

    Args:
        mol: input molecule to protect
        substruct: optional substructure query to identify atoms to protect
        atoms: optional list of atom indices to protect
        in_place: whether to perform the protection in place or return a copy of the molecule
    """
    if atoms is None:
        atoms = []
    elif not isinstance(atoms, (tuple, list)):
        atoms = [atoms]

    # do not perform protection in place
    if in_place:
        mol_copy = mol
    else:
        mol_copy = copy_mol(mol)

    if substruct is not None:
        matches = mol_copy.GetSubstructMatches(substruct)
        atoms.extend(itertools.chain(*matches))

    for a in atoms:
        if a is None:
            continue
        mol_copy.GetAtomWithIdx(a).SetProp("_protected", "1")

    return mol_copy


def add_hs(
    mol: Mol,
    explicit_only: bool = False,
    add_coords: bool = False,
    only_on_atoms: Optional[List[int]] = None,
    add_residue_info: bool = False,
):
    """Adds hydrogens to the molecule.

    Args:
        mol: a molecule.
        explicit_only: whether to only add explicit hydrogens.
        add_coords: whether to add 3D coordinates to the hydrogens.
        only_on_atoms: a list of atoms to add hydrogens only on.
        add_residue_info: whether to add residue information to the hydrogens.
            Useful for PDB files.
    """
    mol = rdmolops.AddHs(
        mol,
        explicitOnly=explicit_only,
        addCoords=add_coords,
        onlyOnAtoms=only_on_atoms,
        addResidueInfo=add_residue_info,
    )

    return mol


def remove_hs(
    mol: Mol,
    implicit_only: bool = False,
    update_explicit_count: bool = False,
    sanitize: bool = True,
):
    """Removes hydrogens from a molecule.

    Args:
        mol: a molecule.
        implicit_only: whether to only remove implicit hydrogens.
        update_explicit_count: whether to update the explicit hydrogen count.
        sanitize: whether to sanitize the molecule after the hydrogens are removed.
    """
    mol = rdmolops.RemoveHs(
        mol,
        implicitOnly=implicit_only,
        updateExplicitCount=update_explicit_count,
        sanitize=sanitize,
    )

    return mol


def strip_mol_to_core(mol: Mol, bond_cutter: Mol = None):
    """Strip a molecule to its core, i.e. remove all atoms not in the core.
    This method 'guess' the molecular core, by finding the ring system.

    Args:
        mol: A molecule.
        bond_cutter: A molecule used to cut the bonds.
    """

    if bond_cutter is None:
        bond_cutter = from_smarts("[R;!$(*=,#[!#6])]!@!=!#[*;$([A;!R][A;!R])]")

    with without_rdkit_log():

        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        out = mol.GetSubstructMatches(bond_cutter)
        bond_inds = [mol.GetBondBetweenAtoms(i, j).GetIdx() for i, j in out]

        if len(bond_inds) > 0:
            fragmented = rdmolops.FragmentOnBonds(mol, bond_inds)
            fragmented = remove_dummies(fragmented)
            fragmented = to_scaffold_murcko(fragmented)
            scaffold = keep_largest_fragment(fragmented)

    return scaffold


def to_scaffold_murcko(mol: Mol, make_generic: bool = False):
    """Extract the Murcko scaffold from a molecule.

    Args:
        mol: A molecule.
        make_generic: Whether to make the scaffold generic.
    """
    scf = MurckoScaffold.GetScaffoldForMol(mol)

    # NOTE(hadim): this is already done in `GetScaffoldForMol`
    # Note sure we need it here.
    scf.UpdatePropertyCache()
    Chem.GetSymmSSSR(scf)  # type: ignore

    if make_generic:
        scf = make_scaffold_generic(scf)
        scf = to_mol(scf)

    return scf


def make_scaffold_generic(mol: Mol, include_bonds: bool = False):
    """Make the atom in a scaffold or molecule generic.

    Args:
        mol: A molecule or a scaffold.
        include_bonds: Whether we should also update bond order or keep as is.
    """

    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 1:
            atom.SetAtomicNum(0)

        atom.SetFormalCharge(0)
        atom.SetChiralTag(rdchem.ChiralType.CHI_UNSPECIFIED)
        atom.SetNoImplicit(0)
        atom.SetNumExplicitHs(0)

    if include_bonds:
        for bond in mol.GetBonds():
            bond.SetBondType(UNSPECIFIED_BOND)

    mol.UpdatePropertyCache()
    Chem.GetSymmSSSR(mol)  # type: ignore

    return mol


def compute_ring_system(mol: Mol, include_spiro: bool = True) -> List[Set[int]]:
    """Compute the list of ring system in a molecule. This is based on RDKit's cookbook:
    https://www.rdkit.org/docs/Cookbook.html#rings-aromaticity-and-kekulization

    Args:
        mol: input molecule
        include_spiro: whether to include spiro rings.

    Returns:
        ring_system: list of ring system (atom indices).
    """
    ri = mol.GetRingInfo()
    systems = []
    for ring in ri.AtomRings():
        ringAts = set(ring)
        nSystems = []
        for system in systems:
            nInCommon = len(ringAts.intersection(system))
            if nInCommon and (include_spiro or nInCommon > 1):
                ringAts = ringAts.union(system)
            else:
                nSystems.append(system)
        nSystems.append(ringAts)
        systems = nSystems
    return systems


def set_atom_positions(
    mol: Mol,
    positions: npt.ArrayLike,
    conf_id: int = 0,
    copy: bool = True,
    use_atom_map_numbers: bool = True,
    remove_previous_conformers: bool = True,
) -> Mol:
    """Add a conformer to a molecule given the atom's positions.

    The conformer 3D flag is automatically set if all the z coordinates are 0.

    **Example:**

    The below example is common when you want to reconstruct a molecule object
    from its SMILES and its held out atomic positions. The position array
    is ordered according to the atom number seen in the SMILES. This is a common
    data structure when working with Quantum Mechanics dataset.

    ```python

    import datamol as dm
    import numpy as np

    # We start with a SMILES where every atoms is mapped to a specific number
    smiles = "[H:14][c:5]1[c:3]([c:7]([c:4]([c:6]([c:8]1[N:10]([H:18])[C:2](=[N+:11]([H:19])[H:20])[N:9]([H:16])[H:17])[H:15])[H:13])[F:1])[H:12]"

    # Every atom position below is mapped to the atom number in the SMILES above
    positions = [
        [1.7, -6.67, 3.15],
        [0.2, 4.72, 0.78],
        [3.54, -2.64, 2.88],
        [0.43, -3.87, -0.09],
        [3.44, -0.2, 1.8],
        [0.02, -1.5, -1.0],
        [2.12, -4.54, 1.9],
        [1.5, 0.48, 0.02],
        [0.53, 7.24, 0.25],
        [1.17, 2.91, -0.85],
        [-1.22, 4.15, 2.71],
        [4.64, -3.24, 4.55],
        [-0.89, -5.43, -0.78],
        [4.52, 1.43, 2.45],
        [-1.45, -1.02, -2.48],
        [-0.15, 8.68, 1.38],
        [1.65, 7.88, -1.21],
        [2.24, 3.64, -2.15],
        [-1.96, 2.4, 3.0],
        [-2.02, 5.59, 3.71],
    ]

    # We build the mol object by setting `remove_hs` to `False`.
    # This is important so the hydrogens and their atom number are preserved.
    mol = dm.to_mol(smiles, remove_hs=False)

    # If you plot the molecule with `dm.to_image(mol)`, you'll notice
    # the atom numbers are added to the drawing.

    # Now we set the atom positions to the newly constructed `mol` object.
    # Here it's important to set `use_atom_map_numbers` to `True`, so the atom numbers
    # from the SMILES are used to match the input positions array.
    # Under the hood, RDKit has set the `molAtomMapNumber` property to all the atoms in the
    # molecule.
    new_mol = dm.set_atom_positions(
        mol=mol,
        positions=positions,
        conf_id=0,
        use_atom_map_numbers=True,
    )

    # The newly set conformer now had the correct 3D positions.
    # You can visualize the molecule with `dm.to_image(new_mol)` in 2D
    # or `dm.viz.conformers(new_mol)`.
    ```

    Args:
        mol: A molecule.
        positions: An array or a list of atomic positions. Shape of `[n_atoms, 3]`.
        conf_id: The conformer ID to set the conformer to.
        copy: Whether to copy the molecule.
        use_atom_map_numbers: Whether to input positions are ordered given the atom mapped
            numbers set in the `molAtomMapNumber` atom property keys. If set to False,
            the default atom indices order is assumed.
        remove_previous_conformers: Whether to remove the previous conformers if any in the input molecule.
    """

    if copy:
        mol = copy_mol(mol)

    # Convert to an array
    positions = np.array(positions)

    # Perform a few sanity checks

    if positions.ndim != 2:
        raise ValueError(f"The positions array is not of dimension 2. Found: {positions.ndim}.")

    excpected_shape = (mol.GetNumAtoms(), 3)
    if positions.shape != excpected_shape:
        raise ValueError(
            f"The shape of `positions` is {positions.shape} instead of {excpected_shape}."
        )

    if use_atom_map_numbers and not all(
        ["molAtomMapNumber" in a.GetPropsAsDict() for a in mol.GetAtoms()]
    ):
        raise ValueError(
            f"The atoms of the input molecule does not contain the molAtomMapNumber property."
            "Set it before calling this function or set `from_atom_map_numbers` to `False`."
        )

    # Remove previous conformers
    if remove_previous_conformers:
        mol.RemoveAllConformers()

    # Remap the rows order in `positions` so it matches
    # with the atom map numbers.
    if use_atom_map_numbers:
        mapped_indices = np.array([a.GetAtomMapNum() for a in mol.GetAtoms()]) - 1
        positions = positions[mapped_indices, :]

    # Check if it's 3D or 2D coords
    is_3d = not np.all(positions[:, 2] == 0)

    # Create the conformer object
    conf = rdchem.Conformer()
    conf.Set3D(is_3d)
    conf.SetId(conf_id)

    # Set the positions
    for i, xyz in enumerate(positions):
        conf.SetAtomPosition(i, rdGeometry.Point3D(*xyz.tolist()))

    # Add the conformer to the molecule
    mol.AddConformer(conf)

    return mol


def get_atom_positions(
    mol: Mol,
    conf_id: int = -1,
    reorder_to_atom_map_number: bool = False,
) -> np.ndarray:
    """Return the atom positions of a given molecule.

    Args:
        mol: A molecule.
        conf_id: The conformer ID to set the conformer to.
        reorder_to_atom_map_number: Whether to reorder the positions to map the
            atom map numbers given by the `molAtomMapNumber` atom property.
    """

    if mol.GetNumConformers() == 0:
        raise ValueError("This molecule does not have conformers.")

    if reorder_to_atom_map_number and not all(
        ["molAtomMapNumber" in a.GetPropsAsDict() for a in mol.GetAtoms()]
    ):
        raise ValueError(
            f"The atoms of the input molecule does not contain the molAtomMapNumber property."
            "Set it before calling this function or set `from_atom_map_numbers` to `False`."
        )

    conformer = mol.GetConformer(id=conf_id)
    positions = conformer.GetPositions()

    if reorder_to_atom_map_number:

        # Remap the rows order in `positions` so it matches
        # with the atom map numbers.
        mapped_indices = np.array([a.GetAtomMapNum() for a in mol.GetAtoms()]) - 1
        positions = positions[mapped_indices, :]

    return positions
