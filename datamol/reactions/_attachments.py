from typing import cast
from typing import Union

import re
import operator

import datamol as dm
from rdkit import Chem

ATTACHMENT_POINT_TOKEN = "*"
ATTACHMENT_POINT_NUM_REGEXP = r"\[{}:?(\d*)\]".format(re.escape(ATTACHMENT_POINT_TOKEN))
ATTACHMENT_POINT_REGEXP = r"(?:{0}|\[{0}[^\]]*\])".format(re.escape(ATTACHMENT_POINT_TOKEN))
ATTACHMENT_POINT_NO_BRACKETS_REGEXP = r"(?<![:\[]){0}(?![:\]])".format(
    re.escape(ATTACHMENT_POINT_TOKEN)
)
ALL_POSSIBLE_ATTACHMENTS = r"\[(\d*){}:?(\d*)\]".format(re.escape(ATTACHMENT_POINT_TOKEN))


def add_brackets_to_attachment_points(smiles: str) -> str:
    """
    Adds brackets to the attachment points (if they don't have them).
    Example: "CC(C)CO*" to "CC(C)CO[*]"

    Args:
        smiles: A smiles string.

    Returns:
        A smiles string with brackets.
    """
    return re.sub(
        ATTACHMENT_POINT_NO_BRACKETS_REGEXP,
        "[{}]".format(ATTACHMENT_POINT_TOKEN),
        smiles,
    )


def convert_attach_to_isotope(
    mol_or_smiles: Union[dm.Mol, str],
    same_isotope: bool = False,
    as_smiles: bool = False,
) -> Union[dm.Mol, str]:
    """Convert attachment to isotope mapping.

    Examples: "O=C(NCc1cnc([*])c1)[*]" to  "O=C(NCc1cnc([1*])c1)[2*]"

    Args:
        mol_or_smiles: A Mol object or a smiles to be converted
        same_isotope: Whether convert to the same isotope.
            Example: "O=C(NCc1cnc([*])c1)[*]" to  "O=C(NCc1cnc([1*])c1)[1*]"

    Returns:
        Converted Mol object or SMILES.
    """
    mol = dm.to_mol(mol_or_smiles)
    smiles = dm.to_smiles(mol)
    smiles = cast(str, smiles)

    smiles = add_brackets_to_attachment_points(smiles)

    # reg matching seems to be the most effective
    subs_reg = r"[\g<1>{}]"
    if same_isotope:
        subs_reg = "[1{}]"

    smiles = re.sub(ATTACHMENT_POINT_NUM_REGEXP, subs_reg.format(ATTACHMENT_POINT_TOKEN), smiles)

    if as_smiles:
        return smiles
    return dm.to_mol(smiles)


def num_attachment_points(mol_or_smiles: Union[dm.Mol, str]) -> int:
    """
    Get the number of attachment point in the

    Args:
        mol_or_smiles: A Mol object or a smiles to be converted

    Returns:
        Number of attachment points of the given molecule.
    """
    if isinstance(mol_or_smiles, dm.Mol):
        mol = cast(dm.Mol, mol_or_smiles)
        n_points = len(
            [atom for atom in mol.GetAtoms() if atom.GetSymbol() == ATTACHMENT_POINT_TOKEN]
        )
    else:
        n_points = len(re.findall(ATTACHMENT_POINT_REGEXP, mol_or_smiles))

    return n_points


def open_attach_points(
    mol: dm.Mol,
    fix_atom_map: bool = False,
    bond_type: dm.BondType = dm.SINGLE_BOND,
) -> dm.Mol:
    """Compute attachment points on a molecule.
    This will highlight all valid attachment point on the current molecule instead.

    Args:
        mol: A Mol object to be processed.
        fix_atom_map: Whether fix the atom mapping of the molecule.
        bond_type: The bond type to be opened.

    Returns:
        Molecule with open attachment points
    """

    emol = Chem.rdchem.RWMol(dm.to_mol(mol))
    with dm.log.without_rdkit_log():
        atoms = [
            (a.GetIdx(), a)
            for a in emol.GetAtoms()
            if a.GetSymbol() != ATTACHMENT_POINT_TOKEN
            and a.GetImplicitValence() > 0
            and (not a.HasProp("_protected") or a.GetProp("_protected") != "1")
        ]
        atoms.sort(reverse=True, key=operator.itemgetter(0))

        for atom in atoms:
            new_atom = Chem.rdchem.Atom(ATTACHMENT_POINT_TOKEN)
            new_atom.SetAtomMapNum(1 if fix_atom_map else atom[0])
            new_index = emol.AddAtom(new_atom)
            emol.UpdatePropertyCache(strict=False)
            if bond_type is not None:
                emol.AddBond(atom[0], new_index, bond_type)
            else:
                emol.AddBond(atom[0], new_index)

    mol = dm.sanitize_mol(emol)
    return mol
