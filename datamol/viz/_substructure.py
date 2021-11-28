import itertools

from typing import Union
from typing import List

import datamol as dm

from ._viz import to_image


def find_bonds_from_atom_list(mol: dm.Mol, atom_indices: List[List[int]]):
    """Return a list of list of existing bond indices between a list
    of list of atom indices.

    Args:
        mol: A molecule.
        atom_indices: A list of list of atom indices.
    """

    bonds = []

    for atom_indices_single in atom_indices:

        bonds_single = []

        # Test all possible combinations for existing bonds
        for a1_idx, a2_idx in itertools.combinations(atom_indices_single, 2):
            bond = mol.GetBondBetweenAtoms(a1_idx, a2_idx)
            if bond is not None:
                bonds_single.append(bond.GetIdx())

        bonds.append(bonds_single)

    return bonds


def match_substructure(
    mols: Union[List[dm.Mol], dm.Mol],
    patterns: Union[List[dm.Mol], dm.Mol],
    highlight_bonds: bool = True,
    **kwargs,
):
    """Generate an image of molecule(s) with substructure matches for a given
    pattern or substructure.

    Args:
        mols: One or more molecules.
        patterns: One or more patterns.
        highlight_bonds: Whether to also highlight the bonds matching the patterns.
        kwargs: Other kwargs passed to `dm.viz.to_image`.
    """

    # NOTE(hadim): `MolsToGridImage` used in `to_image` can't use a list of list of indices
    # for every molecules so it's not really possible to have different colors for different
    # matches in the same molecules.
    # In the future, we will implement our custom `MolsToGridImage` in order to have more controls
    # on the colors used.
    # For the same reason, we don't bother about colors here.

    if isinstance(mols, dm.Mol):
        mols = [mols]

    if isinstance(patterns, dm.Mol):
        patterns = [patterns]

    # Always copy mols and patterns
    mols = [dm.copy_mol(mol) for mol in mols]
    patterns = [dm.copy_mol(mol) for mol in patterns]

    all_atom_indices = []
    all_bond_indices = []

    for mol in mols:

        atom_indices = []
        bond_indices = []

        for pattern in patterns:

            atom_indices += list(mol.GetSubstructMatches(pattern, uniquify=True))

            if highlight_bonds:
                bond_indices += dm.viz.find_bonds_from_atom_list(mol, atom_indices)
            else:
                bond_indices += []

        # NOTE(hadim): we must flatten the atom/bond indices, since `MolsToGridImage`
        # don't accept multiple list of indices for every single molecule.
        bond_indices = [item for sublist in bond_indices for item in sublist]
        atom_indices = [item for sublist in atom_indices for item in sublist]

        all_atom_indices.append(atom_indices)
        all_bond_indices.append(bond_indices)

    image = to_image(
        mols,
        highlight_atom=all_atom_indices,
        highlight_bond=all_bond_indices,
        # highlightAtomColors=atom_colors,
        # highlightBondColors=bond_colors,
        **kwargs,
    )

    return image
