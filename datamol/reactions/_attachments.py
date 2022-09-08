from typing import Union
import re
import datamol as dm
import operator
from rdkit import Chem

ATTACHMENT_POINT_TOKEN = "*"
ATTACHMENT_POINT_NUM_REGEXP = r"\[{}:?(\d*)\]".format(re.escape(ATTACHMENT_POINT_TOKEN))
ATTACHMENT_POINT_REGEXP = r"(?:{0}|\[{0}[^\]]*\])".format(re.escape(ATTACHMENT_POINT_TOKEN))
ATTACHMENT_POINT_NO_BRACKETS_REGEXP = r"(?<![:\[]){0}(?![:\]])".format(
    re.escape(ATTACHMENT_POINT_TOKEN)
)
ALL_POSSIBLE_ATTACHMENTS = r"\[(\d*){}:?(\d*)\]".format(re.escape(ATTACHMENT_POINT_TOKEN))


def add_brackets_to_attachment_points(smi: str):
    """
    Adds brackets to the attachment points (if they don't have them).
    """
    return re.sub(
        ATTACHMENT_POINT_NO_BRACKETS_REGEXP,
        "[{}]".format(ATTACHMENT_POINT_TOKEN),
        smi,
    )


def convert_attach_to_isotope(mol_or_smi: Union[Chem.Mol, str], same_isotope:bool=False):
    """Convert attachment to isotope mapping"""
    mol = dm.to_mol(mol_or_smi)
    smiles = dm.to_smiles(mol)
    smiles = add_brackets_to_attachment_points(smiles)
    # reg matching seems to be the most effective
    subs_reg = "[\g<1>{}]"
    if same_isotope:
        subs_reg = "[1{}]"
    smiles = re.sub(ATTACHMENT_POINT_NUM_REGEXP, subs_reg.format(ATTACHMENT_POINT_TOKEN), smiles)
    return dm.to_mol(smiles)


def num_attachment_points(mol_or_smi):
    """Get the number of attachment point in the"""
    if isinstance(mol_or_smi, Chem.Mol):
        return len(
            [atom for atom in mol_or_smi.GetAtoms() if atom.GetSymbol() == ATTACHMENT_POINT_TOKEN]
        )
    return len(re.findall(ATTACHMENT_POINT_REGEXP, mol_or_smi))


def open_attach_points(mol, fix_atom_map=False, bond_type=dm.SINGLE_BOND):
    """Compute attachment points on a molecule
    This will highlight all valid attachment point on the current molecule instead
    """
    emol = Chem.RWMol(dm.to_mol(mol))
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
            new_atom = Chem.Atom(ATTACHMENT_POINT_TOKEN)
            new_atom.SetAtomMapNum(1 if fix_atom_map else atom[0])
            new_index = emol.AddAtom(new_atom)
            emol.UpdatePropertyCache(strict=False)
            if bond_type is not None:
                emol.AddBond(atom[0], new_index, bond_type)
            else:
                emol.AddBond(atom[0], new_index)
    mol = dm.sanitize_mol(emol)
    return mol
