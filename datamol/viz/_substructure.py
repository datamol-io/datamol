from typing import Any
from typing import Union
from typing import List

import datamol as dm

from ._viz import to_image


def match_substructure(
    mols: Union[List[dm.Mol], dm.Mol],
    queries: Union[List[dm.Mol], dm.Mol],
    highlight_bonds: bool = True,
    copy: bool = True,
    **kwargs: Any,
):
    """Generate an image of molecule(s) with substructure matches for a given
    pattern or substructure.

    Args:
        mols: One or more molecules.
        queries: One or more queries.
        highlight_bonds: Whether to also highlight the bonds matching the patterns.
        copy: Whether to copy the molecules and the queries.
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

    if isinstance(queries, dm.Mol):
        queries = [queries]

    # Copy mols and patterns
    if copy:
        mols = [dm.copy_mol(mol) for mol in mols]
        queries = [dm.copy_mol(mol) for mol in queries]

    all_atom_indices = []
    all_bond_indices = []

    for mol in mols:

        atom_indices = []
        bond_indices = []

        for query in queries:
            if highlight_bonds:
                atom_matches, bond_matches = dm.substructure_matching_bonds(mol, query)
                atom_indices += atom_matches
                bond_indices += bond_matches
            else:
                atom_indices += list(mol.GetSubstructMatches(query, uniquify=True))  # type: ignore
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
        **kwargs,
    )

    return image
