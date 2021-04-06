import copy
import itertools
import operator
import random

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolops

import datamol as dm


def pick_atom_idx(mol, prepick=None):
    """pick an atom from the molecule"""
    mol.UpdatePropertyCache()
    if not (prepick is not None and prepick >= 0 and prepick < mol.GetNumAtoms()):
        pickable_atoms = [a.GetIdx() for a in mol.GetAtoms() if a.GetImplicitValence() > 0]
        if pickable_atoms:
            prepick = random.choice(pickable_atoms)
        else:
            prepick = None
    return prepick


def add_bond_between(mol, a1, a2, bond_type):
    """Add a new bond between atom"""
    emol = Chem.EditableMol(mol)
    emol.AddBond(a1.GetIdx(), a2.GetIdx(), bond_type)
    return dm.sanitize_mol(emol.GetMol())


def update_bond(mol, bond, bond_type):
    """Update bond type between atoms"""
    new_mol = dm.copy_mol(mol)
    with dm.without_rdkit_log():
        new_bond = new_mol.GetBondWithIdx(bond.GetIdx())
        new_bond.SetBondType(bond_type)
    return dm.sanitize_mol(new_mol)


def _all_atom_join(mol, a1, a2):
    """Join two atoms (a1, a2) in a molecule in all possible valid manner"""
    new_mols = []
    with dm.without_rdkit_log():
        try:
            Chem.Kekulize(mol, clearAromaticFlags=True)
        except:
            pass
        v1, v2 = a1.GetImplicitValence(), a2.GetImplicitValence()
        bond = mol.GetBondBetweenAtoms(a1.GetIdx(), a2.GetIdx())
        if bond is None:
            if v1 > 0 and v2 > 0:
                new_mols.append(add_bond_between(mol, a1, a2, dm.SINGLE_BOND))
            if v1 > 1 and v2 > 1:
                new_mols.append(add_bond_between(mol, a1, a2, dm.DOUBLE_BOND))
            if v1 > 2 and v2 > 2:
                new_mols.append(add_bond_between(mol, a1, a2, dm.TRIPLE_BOND))

        elif bond.GetBondType() == dm.SINGLE_BOND:
            if v1 > 0 and v2 > 0:
                new_mols.append(update_bond(mol, bond, dm.DOUBLE_BOND))
            if v1 > 1 and v2 > 1:
                new_mols.append(update_bond(mol, bond, dm.TRIPLE_BOND))

        elif bond.GetBondType() == dm.DOUBLE_BOND:
            if v1 > 0 and v2 > 0:
                new_mols.append(update_bond(mol, bond, dm.TRIPLE_BOND))
    return [mol for mol in new_mols if mol is not None]


def _compute_fragment_join(
    mol,
    fragment,
    mol_atom_count,
    bond_between_rings=True,
    asMols=True,
):
    """List all posibilities of where a fragment can be attached to a mol"""
    fragment = copy.copy(
        fragment
    )  # need to copy the fragment copy is faster than all the other methods
    with dm.without_rdkit_log():
        combined = Chem.CombineMols(mol, fragment)
        for i1 in range(mol.GetNumAtoms()):
            a1 = combined.GetAtomWithIdx(i1)
            if a1.GetImplicitValence() == 0:
                continue
            for i2 in range(fragment.GetNumAtoms()):
                i2 += mol_atom_count
                a2 = combined.GetAtomWithIdx(i2)
                if a2.GetImplicitValence() == 0:
                    continue
                # no bond between atoms already in rings
                if not bond_between_rings and a1.IsInRing() and a2.IsInRing():
                    continue
                # no bond to form large rings
                else:
                    possibilities = _all_atom_join(combined, a1, a2)
                    for x in possibilities:
                        x = dm.sanitize_mol(x)
                        if x is not None:
                            if not asMols:
                                x = dm.to_smiles(x)
                            yield x


def _compute_mmpa_assembly(cores, side_chains, max_num_action=float("Inf")):
    """Enumerate core and side_chains assembly combination.
    Input Core and side_chain are expected to have [1*] in place of the attachment point

    Note that this is based on a dm.SINGLE_BOND mmpa cut.

    Arguments
    ----------
        cores: list of <Chem.Mol>
            List of core
        side_chains: list of <Chem.Mol>
            List of side chains
        max_num_action: int, optional
            Maximum number of assembly
            (Default: inf)

    Returns
    -------
        res: list of <Chem.Mol>
            Molecules obtained by merging core and side_chains
    """
    reaction = AllChem.ReactionFromSmarts("[*:1]-[1*].[1*]-[*:2]>>[*:1]-[*:2]")
    molecules = []
    n_seen = 0
    random.shuffle(side_chains)
    stop = False
    for core in cores:
        if not stop:
            for sidechain in side_chains:
                molecules.append(reaction.RunReactants((core, sidechain))[0][0])  # first only
                n_seen += 1
                if n_seen > max_num_action:
                    stop = True
                    break
    random.shuffle(molecules)
    return molecules


def all_join_on_attach_point(mol1, mol2):
    """Join two molecules on all possible attaching point

    Arguments
    ---------
        mol1: <Chem.Mol>
            input molecule 1
        mol2: <Chem.Mol>
            input molecule 2

    Returns
    -------
        iterator of all possible way to attach both molecules from dummy indicators.
    """
    atom_map_min = 100
    mol_idxs = []
    count = 0
    mod_mols = []
    for ind, m in enumerate([mol1, mol2]):
        atms = [(a.GetIdx(), a) for a in m.GetAtoms() if not a.IsInRing() and a.GetAtomicNum() == 0]
        atms.sort(reverse=True, key=operator.itemgetter(0))
        for a_idx, a in atms:
            for a_nei in a.GetNeighbors():
                a_nei.SetAtomMapNum(atom_map_min + count)
                count += 1
        mod_mol = dm.fix_mol(m)
        mod_mols.append(mod_mol)
        mol_idxs.append(
            [a.GetIdx() for a in mod_mol.GetAtoms() if a.GetAtomMapNum() >= atom_map_min]
        )
    for ind1, ind2 in itertools.product(*mol_idxs):
        yield random_fragment_add(copy.copy(mod_mols[0]), copy.copy(mod_mols[1]), ind1, ind2)


def all_fragment_attach(
    mol,
    fragmentlist,
    bond_between_rings=True,
    max_num_action=10,
    asMols=True,
):
    """List all possible way to attach a list of fragment to a dm.SINGLE_BOND molecule.

    .. warning::
        This is computationally expensive

    Arguments
    ----------
        mol: <Chem.Mol>
            Input molecule
        fragmentlist: list of <Chem.Mol>
            Molecular fragments to attach
        bond_between_rings: bool, optional
            Whether to allow bond between two rings atoms
            (Default: True)
        max_num_action: int, optional
            Maximum fragment attachment to allow. Reduce time complexity
            (Default: 10)
        asMols: bool, optional
            Whether to return output as molecule or smiles
    Returns
    -------
        All possible molecules resulting from attaching the molecular fragment to the root molecule

    """
    fragment_set = set([])
    mol_atom_count = mol.GetNumAtoms()
    generators = [None] * len(fragmentlist)
    empty_generators = np.zeros(len(generators))
    while len(fragment_set) < max_num_action and not np.all(empty_generators):
        for i, fragment in enumerate(fragmentlist):
            if len(fragment_set) >= max_num_action:
                break
            if generators[i] is None:
                generators[i] = _compute_fragment_join(
                    mol, fragment, mol_atom_count, bond_between_rings, asMols
                )
            if not empty_generators[i]:
                try:
                    fragment_set.add(next(generators[i]))
                except StopIteration as e:
                    empty_generators[i] = 1
                    continue
    return fragment_set


def all_atom_add(
    mol,
    atom_types=["C", "N", "O", "F", "Cl", "Br"],
    asMols=True,
    max_num_action=float("Inf"),
    **kwargs,
):
    """Add a new atom on the mol, by considering all bond type

    .. warning::
        This is computationally expensive

    Arguments
    ----------
        mol: <Chem.Mol>
            Input molecule
        atom_types: list
            List of atom symbol to use as replacement
            (Default: ["C", "N", "O", "F", "Cl", "Br"])
        asMols: bool, optional
            Whether to return output as molecule or smiles
        max_num_action: float, optional
            Maximum number of action to reduce complexity
    Returns
    -------
        All possible molecules with one additional atom added

    """
    new_mols = []
    stop = False
    with dm.without_rdkit_log():
        for atom in mol.GetAtoms():
            if stop:
                break
            if atom.GetImplicitValence() == 0:
                continue
            for atom_symb in atom_types:
                emol = Chem.RWMol(mol)
                new_index = emol.AddAtom(Chem.Atom(atom_symb))
                emol.UpdatePropertyCache(strict=False)
                new_mols.extend(_all_atom_join(emol, atom, emol.GetMol().GetAtomWithIdx(new_index)))
                if len(new_mols) > max_num_action:
                    stop = True
                    break

        new_mols = [dm.sanitize_mol(mol) for mol in new_mols]
        new_mols = [mol for mol in new_mols if mol is not None]
        if not asMols:
            return [dm.to_smiles(x) for x in new_mols if x]
    return new_mols


def all_atom_replace(
    mol, atom_types=["C", "N", "S", "O"], asMols=True, max_num_action=float("Inf"), **kwargs
):
    """Replace all non-hydrogen atoms by other possibilities.

    .. warning::
        This is computationally expensive

    Arguments
    ----------
        mol: <Chem.Mol>
            Input molecule
        atom_types: list
            List of atom symbol to use as replacement
            (Default: ['C', 'N', 'S', 'O'])
        asMols: bool, optional
            Whether to return output as molecule or smiles
        max_num_action: float, optional
            Maximum number of action to reduce complexity

    Returns
    -------
        All possible molecules with atoms replaced

    """
    new_mols = []
    stop = False
    with dm.without_rdkit_log():
        for atom in mol.GetAtoms():
            if stop:
                break
            if atom.GetAtomicNum() > 1:
                for atom_symb in atom_types:
                    emol = Chem.RWMol(mol)
                    emol.ReplaceAtom(atom.GetIdx(), Chem.Atom(atom_symb))
                    new_mols.append(emol)
                    if len(new_mols) > max_num_action:
                        stop = True
                        break

        # Sanitize and remove bad molecules
        new_mols = [dm.sanitize_mol(mol) for mol in new_mols]
        new_mols = [mol for mol in new_mols if mol is not None]

    if not asMols:  # Return SMILES
        return [dm.to_smiles(x) for x in new_mols]
    return new_mols


def all_bond_add(
    mol,
    allowed_ring_sizes=None,
    bond_between_rings=True,
    asMols=True,
    max_num_action=float("Inf"),
    **kwargs,
):
    """Add bond to a molecule

    .. warning::
        This is computationally expensive

    Arguments
    ----------
        mol: <Chem.Mol>
            Input molecule
        allowed_ring_sizes: list, optional
            Set of integer allowed ring sizes; used to remove some
            actions that would create rings with disallowed sizes.
        bond_between_rings: bool, optional
            Whether to allow actions that add bonds
            between atoms that are both in rings.
        asMols: bool, optional
            Whether to return output as molecule or smiles
        max_num_action: float, optional
            Maximum number of action to reduce complexity

    Returns
    -------
        All possible molecules with additional bond added between atoms
    """
    new_mols = []
    num_atoms = mol.GetNumAtoms()
    stop = False
    for i1 in range(num_atoms):
        if stop:
            break
        a1 = mol.GetAtomWithIdx(i1)
        if a1.GetImplicitValence() == 0:
            continue
        for i2 in range(i1 + 1, num_atoms):
            a2 = mol.GetAtomWithIdx(i2)
            # Chem.rdmolops.GetShortestPath(mol, i1, i2)
            all_paths = get_all_path_between(mol, i1, i2, ignore_cycle_basis=True)
            all_path_len = {len(path) for path in all_paths}
            if a2.GetImplicitValence() == 0:
                continue
            # no bond between atoms already in rings
            bond = mol.GetBondBetweenAtoms(i1, i2)
            if not bond_between_rings and a1.IsInRing() and a2.IsInRing():
                continue
            # no bond to form large rings
            if (
                (bond is None)
                and (allowed_ring_sizes is not None)
                and not all_path_len.issubset(allowed_ring_sizes)
            ):
                continue
            new_mols.extend(_all_atom_join(mol, a1, a2))
            if len(new_mols) > max_num_action:
                stop = True
                break
    if not asMols:
        return list({dm.to_smiles(x) for x in new_mols if x})
    return [m for m in new_mols if m is not None]


def all_bond_remove(
    mol: Chem.rdchem.Mol,
    as_mol: bool = True,
    allow_bond_decrease: bool = True,
    allow_atom_trim: bool = True,
    max_num_action=float("Inf"),
):
    """Remove bonds from a molecule

    Warning:
        This can be computationally expensive.

    Args:
        mol: Input molecule
        allow_bond_decrease: Allow decreasing bond type in addition to bond cut
        max_num_action: Maximum number of action to reduce complexity
        allow_atom_trim: Allow bond removal even when it results in dm.SINGLE_BOND

    Returns:
        All possible molecules from removing bonds

    """
    new_mols = []

    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except:
        pass

    for bond in mol.GetBonds():
        if len(new_mols) > max_num_action:
            break

        original_bond_type = bond.GetBondType()
        emol = Chem.RWMol(mol)
        emol.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        new_mol = dm.sanitize_mol(emol.GetMol())

        if not new_mol:
            continue

        frag_list = list(rdmolops.GetMolFrags(new_mol, asMols=True))
        has_single_atom = any([x.GetNumAtoms() < 2 for x in frag_list])
        if not has_single_atom or allow_atom_trim:
            new_mols.extend(frag_list)
        if allow_bond_decrease:
            if original_bond_type in [dm.DOUBLE_BOND, dm.TRIPLE_BOND]:
                new_mol = update_bond(mol, bond, dm.SINGLE_BOND)
                if new_mol is not None:
                    new_mols.extend(list(rdmolops.GetMolFrags(new_mol, asMols=True)))
            if original_bond_type == dm.TRIPLE_BOND:
                new_mol = update_bond(mol, bond, dm.DOUBLE_BOND)
                if new_mol is not None:
                    new_mols.extend(list(rdmolops.GetMolFrags(new_mol, asMols=True)))

    new_mols = [mol for mol in new_mols if mol is not None]

    if not as_mol:
        return [dm.to_smiles(x) for x in new_mols if x]

    return new_mols


def all_fragment_on_bond(mol, asMols=False, max_num_action=float("Inf"), break_aromatic=True):
    """Fragment all possible bond in a molecule and return the set of resulting fragments
    This is similar to `random_bond_cut`, but is not stochastic as it does not return a random fragment
    but all the fragments resulting from all potential bond break in the molecule.

    .. note::
        This will always be a subset of all_bond_remove, the main difference being that all_bond_remove, allow decreasing
        bond count, while this one will always break a molecule into two.

    Arguments
    ----------
        mol: <Chem.Mol>
            input molecule
        asMols: bool, optional
            Whether to return results as mols or smiles
        max_num_action: float, optional
            Maximum number of action to reduce complexity
        break_aromatic: bool, optional
            Whether to attempt to break even aromatic bonds
            (Default: True)

    Returns
    -------
        set of fragments

    """
    mol.GetRingInfo().AtomRings()
    fragment_set = set([])
    bonds = list(mol.GetBonds())
    stop = False
    if bonds:
        if break_aromatic:
            Chem.Kekulize(mol, clearAromaticFlags=True)
        for bond in bonds:
            if stop:
                break
            if break_aromatic or not bond.GetIsAromatic():
                truncate = Chem.FragmentOnBonds(mol, [bond.GetIdx()], addDummies=False)
                truncate = dm.sanitize_mol(truncate)
                if truncate is not None:
                    for frag in rdmolops.GetMolFrags(truncate, asMols=True):
                        frag = dm.sanitize_mol(frag)
                        if frag:
                            if not asMols:
                                frag = dm.to_smiles(frag)
                            fragment_set.add(frag)
                        if len(fragment_set) > max_num_action:
                            stop = True
                            break
    return fragment_set


def all_fragment_update(
    molparent,
    fragmentlist,
    bond_between_rings=True,
    max_num_action=float("Inf"),
    asMols=False,
    **kwargs,
):
    """
    Break molecule a molecules into all set of fragment (including the molecule itself).
    Then enumerate all possible combination with blocks from the fragmentlist.
    This corresponds to exploring all valid actions by adding/replacing fragments in a molecules.

    .. warning::
        This is computationally expensive

    .. note::
        You should perform a valency check after

    Arguments
    ----------
        molparent: <Chem.Mol>
            input molecule
        fragmentlist: list
            List of blocks to use for replacement, or addition to molparent
        bond_between_rings: bool, optional
            Whether to allow bond between rings
            (Default: True)
        max_num_action: float, optional
            Maximum number of action to reduce complexity
        asMols: bool, optional
            Whether to return smiles or mols

    Returns
    -------
        set of modified mols
    """
    fragment_set = set([])
    mol_frags = anybreak(molparent, rem_parent=False)
    for mol in mol_frags:
        mol_update = all_fragment_attach(
            mol, fragmentlist, bond_between_rings, max_num_action, asMols
        )
        fragment_set.update(mol_update)
        if len(fragment_set) > max_num_action:
            break
    return list(fragment_set)


def all_mmpa_assemble(molist, max_num_action=float("Inf"), asMols=True, **kwargs):
    """Enumerate all mmpa assembly of molecules in molist

    Arguments
    ----------
        molist: list of <Chem.Mol>
            List of molecules to fragmente and reconstruct
        asMols: bool, optional
            Whether to return smiles or mols
        max_num_action: int, optional
            Maximum number of assembly
            (Default: inf)

    Returns
    -------
        res: list of <Chem.Mol>
            Molecules obtained by merging core and side_chains
    """
    frags = set([])
    cores = []
    side_chains = []
    for mol in molist:
        mol_frag = mmpa_frag(mol, max_bond_cut=30)
        if not mol_frag:
            continue
        _, mol_frag = map(list, zip(*mol_frag))
        for m in mol_frag:
            core, sidechain = m.split(".")
            cores.append(Chem.MolFromSmiles(core.replace("[*:1]", "[1*]")))
            side_chains.append(Chem.MolFromSmiles(sidechain.replace("[*:1]", "[1*]")))
    new_mols = _compute_mmpa_assembly(cores, side_chains, max_num_action=max_num_action)
    if not asMols:
        new_mols = [dm.to_smiles(x) for x in new_mols if x]
    return new_mols


def all_fragment_assemble(
    fragmentlist,
    max_num_action=float("Inf"),
    asMols=True,
    seen=None,
    **kwargs,
):
    """Assemble a set of fragment into a new molecule

    .. warning::
        This is computationally expensive

    Arguments
    ----------
        fragmentlist: list
            List of blocks to use for replacement, or addition to molparent
        max_num_action: float, optional
            Maximum number of action to reduce complexity. No limit by default
        asMols: bool, optional
            Whether to return smiles or mols
        seen: list, optional
            List of initial molecules

    Returns
    -------
        reconstructed molecules

    """
    mols = []
    for m in dm.assemble.assemble_brics_order(
        fragmentlist, seen=seen, allow_incomplete=False, max_n_mols=max_num_action
    ):
        if len(mols) > max_num_action:
            break
        mols.append(m)

    if not asMols:
        mols = [dm.to_smiles(x) for x in mols if x is not None]
    return mols


def all_transform_apply(
    mol,
    rxns,
    max_num_action=float("Inf"),
    asMols=True,
    **kwargs,
):
    """
    Apply a transformation defined as a reaction from a set of reaction to the input molecule.

    The reaction need to be one reactant-only

    Arguments
    ----------
        mol: <Chem.Mol>
            Input molecule
        rnxs: list
            list of reactions/ reaction smarts
        max_num_action: int, optional
            Maximum number of result to return
            (Default: inf)
        asMols: bool, optional
            Whether to return smiles or mols

    Returns
    -------
        Products obtained from applying the chemical reactions
    """

    mols = set([])
    with dm.without_rdkit_log():
        for rxn in rxns:
            if len(mols) >= max_num_action:
                break
            if isinstance(rxn, str):
                rxn = AllChem.ReactionFromSmarts(rxn)
            try:
                pcdts = [products[0] for products in rxn.RunReactants([mol])]
                pcdts = [dm.sanitize_mol(x) for x in pcdts]
                mols.update([dm.to_smiles(x) for x in pcdts if x])
            except:
                pass
    mols = [x for x in mols if x is not None]
    if np.isfinite(max_num_action):
        mols = mols[:max_num_action]

    mols = [dm.to_mol(x) for x in mols]
    if not asMols:
        mols = [dm.to_smiles(x) for x in mols if x is not None]
    return mols


def mmpa_fragment_exchange(mol1, mol2, return_all=False, **kwargs):
    """Perform a fragment exchange between two molecules using mmpa rules

    Arguments
    ----------
        mol1: <Chem.Mol>
            input molecule 1
        mol2: <Chem.Mol>
            input molecule 1
        return_all: bool, optional
            Whether to return list of all molecules

    Returns
    -------
        modified_mol1, modified_mol2
            Molecules obtained by exchanging fragment between mol1 and mol2.
            In case of failure, mol1, mol2 are returned

    """

    unwanted = [dm.to_smiles(m) for m in [mol1, mol2]] + [None]
    res = all_mmpa_assemble([mol1, mol2])
    # find unique
    res = set([dm.to_smiles(m) for m in res])
    res = list(res - set(unwanted))
    out = []
    for sm in res:
        r = None
        try:
            r = dm.to_mol(sm, sanitize=True)
        except:
            continue
        if r is not None:
            out.append(r)

    if return_all:
        return out
    random.shuffle(out)
    out.extend([mol1, mol2])
    return out[0], out[1]
