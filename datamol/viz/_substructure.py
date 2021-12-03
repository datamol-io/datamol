import itertools

from typing import Union
from typing import List

import datamol as dm

from ._viz import to_image


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

            matches = list(mol.GetSubstructMatches(pattern, uniquify=True))

            if highlight_bonds:
                bond_indices += [
                    dm.atom_list_to_bond(mol, match, bond_as_idx=True) for match in matches
                ]
            else:
                bond_indices += []

            atom_indices += matches

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
