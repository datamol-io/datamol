import copy
import itertools
import operator
import random

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolops

import datamol as dm


def _pick_atom_idx(mol, prepick=None):
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
    new_mol = copy.copy(mol)
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
    mol,
    asMols=True,
    allow_bond_decrease=True,
    allow_atom_trim=True,
    max_num_action=float("Inf"),
):
    """
    Remove bonds from a molecule

    .. warning::
        This is could be computationally expensive

    Arguments
    ----------
        mol: <Chem.Mol>
            Input molecule
        allow_bond_decrease: bool, optional
            Allow decreasing bond type in addition to bond cut
            (Default: True)
        max_num_action: float, optional
            Maximum number of action to reduce complexity
        allow_atom_trim: bool, optional
            Allow bond removal even when it results in dm.SINGLE_BONDton atoms
            (Default: True)


    Returns
    -------
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
        has_dm.SINGLE_BOND_atom = any([x.GetNumAtoms() < 2 for x in frag_list])
        if not has_dm.SINGLE_BOND_atom or allow_atom_trim:
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
    if not asMols:
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


def random_transform_apply(mol, rxns, probs=None, **kwargs):
    """Apply a random transform on a molecule and return product

    Arguments
    ----------
        mol: <Chem.Mol>
            Input molecule
        rnxs: list
            list of reactions/ reaction smarts
        probs: list, optional
            probability of applying each operation in rxns
            (Default: None)

    Returns
    --------
        mutated_mol: <Chem.Mol>
            randomly mutated molecule
    """
    rxn = random.choices(rxns, k=1, weights=probs)
    res = all_transform_apply(mol, rxn)
    with dm.without_rdkit_log():
        res = dm.sanitize_best(res) or mol
    return res


def random_bond_decrease(
    mol,
    bond=None,
    largest_fragment=True,
    kekulize=False,
    **kwargs,
):
    """
    Randomly change bond types (dm.SINGLE_BOND, dm.DOUBLE_BOND or dm.TRIPLE_BOND) by decreasing the bond count between a given pair of atoms
    For aromatic bonds, it's a bit tricky and kekulization is attempted first.

    Arguments
    ----------
        mol: <Chem.Mol>
            Input molecule
        bond: list of <Chem.Bond>
            Use this specifc bond from the input molecule if provided. Otherwise, randomly select a bond.
            (Default: None)
        largest_fragment: bool, optional
            Return largest fragment is set to True if molecule was broken into fragments
            (Default: True)
        kekulize: bool, optional
            Attempt molecule kekulization first to ensure aromatic bond are not present.
            (Default: True)

    Returns
    -------
        Return a molecule with a random bond modification or original mol when this fails
    """
    if kekulize:
        mol = copy.copy(mol)
        Chem.Kekulize(mol, clearAromaticFlags=False)
    if bond is not None:
        bonds = [bond]
    else:
        bonds = [
            bond
            for bond in mol.GetBonds()
            if bond.GetBondType() in [dm.SINGLE_BOND, dm.DOUBLE_BOND, dm.TRIPLE_BOND]
        ]
    if bonds:
        bond = random.choice(bonds)
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        bondtype = dm.decrease_bond(bond)
        if bondtype is not None:
            mol = update_bond(mol, bond, bondtype) or mol
        else:
            em = Chem.EditableMol(mol)
            em.RemoveBond(a1.GetIdx(), a2.GetIdx())
            mol = dm.sanitize_mol(em.GetMol()) or mol
    if largest_fragment:
        return dm.keep_largest_fragment(mol)
    return mol


def random_atom_add(
    mol,
    atom_types=["C", "N", "O", "F", "Cl", "Br"],
    largest_fragment=True,
    **kwargs,
):
    """Randomly add a new atom on the mol, by considering all bond type.

    Arguments
    ----------
        mol: <Chem.Mol>
            Input molecule
        atom_types: list
            List of atom symbol to use as replacement
            (Default: ["C", "N", "O", "F", "Cl", "Br"])
        largest_fragment: bool, optional
            Return largest fragment is set to True if molecule was broken into fragments
            (Default: True)
    Returns
    -------
        A new molecule with a new atom randomly added or original molecule if failed

    """
    with dm.without_rdkit_log():
        atoms = [a for a in mol.GetAtoms() if a.GetImplicitValence() > 0]
        random.shuffle(atoms)
        i = 0
        stop = False
        while (not stop) and (i < len(atoms)):
            atom = atoms[i]
            i += 1
            random.shuffle(atom_types)
            for atom_symb in atom_types:
                emol = Chem.RWMol(mol)
                new_index = emol.AddAtom(Chem.Atom(atom_symb))
                emol.UpdatePropertyCache(strict=False)
                new_mols = list(_all_atom_join(emol, atom, emol.GetMol().GetAtomWithIdx(new_index)))
                random.shuffle(new_mols)
                out = dm.sanitize_best(new_mols)
                if out is not None:
                    mol = out
                    stop = True
                    break
    if largest_fragment:
        return dm.keep_largest_fragment(mol)
    return mol


def random_atom_replace(
    mol,
    atom_types=["C", "N", "S", "O"],
    largest_fragment=True,
    **kwargs,
):
    """Randomly replace a non-hydrogen atoms by other possibilities.

    .. warning::
        This is computationally expensive

    Arguments
    ----------
        mol: <Chem.Mol>
            Input molecule
        atom_types: list
            List of atom symbol to use as replacement
            (Default: ['C', 'N', 'S', 'O'])
        largest_fragment: bool, optional
            Return largest fragment is set to True if molecule was broken into fragments
            (Default: True)

    Returns
    -------
        Molecule with atom randomly replaced or original molecule

    """

    with dm.without_rdkit_log():
        atoms = [a for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
        random.shuffle(atoms)
        i = 0
        stop = False
        while (not stop) and (i < len(atoms)):
            atom = atoms[i]
            i += 1
            random.shuffle(atom_types)
            for atom_symb in atom_types:
                emol = Chem.RWMol(mol)
                emol.ReplaceAtom(atom.GetIdx(), Chem.Atom(atom_symb))
                out = dm.sanitize_mol(emol)
                if out is not None:
                    mol = out
                    stop = True
                    break
    if largest_fragment:
        return dm.keep_largest_fragment(mol)
    return mol


def random_fragment_add(
    mol,
    fragment,
    mol_atom_idx=None,
    frag_atom_idx=None,
    largest_fragment=True,
    **kwargs,
):
    """
    Given two rdkit mol objects, combine them by adding a bond between two atoms.
    If one of the index is not provided, it will be randomly selected from the set of possible atoms that can be connected
    Note that a sanity check should be performed later. In the worst case, the original molecule will be returned

    Arguments
    ----------
        mol: RDKit Mol
            original molecule
        fragment: RDKit Mol
            Fragment to add to mol
        mol_atom_idx: int, optional
            Atom with the source bond in the molecule
        frag_atom_idx: int, optional
            Second atom in the fragment that should be linked
        largest_fragment: bool, optional
            Whether to keep only largest fragment in the molecule

    Returns
    -------
        A modified molecules or the input molecule if failed
    """
    fragment = copy.copy(fragment)  # need to always copy the fragment
    combined = Chem.CombineMols(mol, fragment)
    if mol.GetNumAtoms() == 0:
        return dm.sanitize_mol(fragment)
    mol_atom_idx = _pick_atom_idx(mol, mol_atom_idx)
    frag_atom_idx = _pick_atom_idx(fragment, frag_atom_idx)
    if mol_atom_idx is not None and frag_atom_idx is not None:
        frag_atom_idx = frag_atom_idx + mol.GetNumAtoms()  # re-adjust the fragment index
        rw_combined = _all_atom_join(
            combined,
            combined.GetAtomWithIdx(mol_atom_idx),
            combined.GetAtomWithIdx(frag_atom_idx),
        )
        rw_combined = [x for x in rw_combined if x is not None]
        random.shuffle(rw_combined)
        rw_combined = dm.sanitize_best(rw_combined + [mol])
        if largest_fragment:
            return dm.keep_largest_fragment(rw_combined)
        return rw_combined
    return mol  # Do nothing !


def random_fragment_del(mol, largest_fragment=True, **kwargs):
    """Remove a fragment from the orignal molecule by
        - first focusing on bond that could breaks the molecule into two fragments
        - then BRICS bond (synthetically available) bond when the former fails

    .. note::
        The probability or returning a fragment is proportional to the size of the fragment.

    Arguments
    ----------
        mol: <Chem.Mol>
            input molecule
        largest_fragment: bool, optional
            Whether to keep only the largest fragment in case molecule is still fragmented
            (Default: True)

    Returns
    -------
        modified molecule or original molecule if failed
    """
    frags = dm.fragment.anybreak(mol)
    frags = [x for x in frags if x is not None]
    if len(frags) > 0:
        truncate = random.choices(frags, k=1, weights=[f.GetNumAtoms() for f in frags])
        mol = dm.sanitize_best(truncate) or mol
    if largest_fragment:
        return dm.keep_largest_fragment(mol)
    return mol


def random_bond_cut(mol, largest_fragment=True, **kwargs):
    """Randomly cut a bond in the input molecule and return the resulting molecule.
    This could either result in fragment deletion, acylisation or molecular fragmentation.

    .. note::
        This could potentially result in molecule for which sanitization fail, since there is
        almost no efficient way of getting away with cutting an aromatic bond.

    Arguments
    ----------
        mol: <Chem.Mol>
            input molecule
        largest_fragment: bool, optional
            Whether to keep only the largest fragment or return full molecule
            (Default: True)

    Returns
    -------
        modified molecule or original molecule if failed
    """
    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except:
        pass
    mol.GetRingInfo().AtomRings()
    bonds = list(mol.GetBonds())
    random.shuffle(bonds)
    if bonds:
        new_mol = None
        while new_mol is None and len(bonds) > 0:
            bond = bonds.pop(0)
            new_mol = Chem.FragmentOnBonds(mol, [bond.GetIdx()], addDummies=False)
            new_mol = dm.sanitize_mol(new_mol)
        mol = new_mol or mol
    if largest_fragment:
        return dm.keep_largest_fragment(mol)
    return mol


def random_fragment_replace(
    mol,
    fragment,
    retry=3,
    largest_fragment=True,
    **kwargs,
):
    """Select a random fragment of a molecule and replace it by a new fragment.
    Replacement does not follow any particular rules except basic chemistry.
    If you want to use attachment points, please see `moms.utils.mols.join_on_attach_point`.

    Performing  this action is tricky in general. Here we use a simple heuristic by breaking the molecule only
    on accessible bond then replacing one of the fragment by a new one.

    Arguments
    ----------
        mol: <Chem.Mol>
            input molecule
        fragment: <Chem.Mol>
            fragment to use to replace part of mol
        retry: int, optional
            number of time to retry until a result is obtained
            (Default: 3)
        largest_fragment: bool, optional
            Whether to keep only the largest fragment or return full molecule
            (Default: True)

    Returns
    -------
        modified molecule or original molecule if failed
    """
    frags = dm.fragment.anybreak(mol)
    if frags:
        new_mol = None
        while retry > 0 and new_mol is None:
            try:
                random.shuffle(frags)
                rand_frag = frags.pop(0)
                new_mols = [
                    dm.sanitize_mol(tmp_mol)
                    for tmp_mol in all_fragment_attach(rand_frag, [fragment], max_num_action=10)
                ]
                new_mols = [tmp_mol for tmp_mol in new_mols if tmp_mol]
                if len(new_mols) > 0:
                    new_mol = random.choice(new_mols)
                    retry = 0
                else:
                    raise ValueError("Frag not valid, next")
            except Exception as e:
                retry -= 1
        mol = new_mol or mol
    if largest_fragment:
        return dm.keep_largest_fragment(mol)
    return mol


def random_fragment_mutation(mol, fragmentlist, largest_fragment=True, **kwargs):
    """Mol mutation performing either fragment replacement or addition

    Arguments
    ----------
        mol: <Chem.Mol>
            input molecule
        fragment: <Chem.Mol>
            fragment to use to replace part of mol
        retry: int, optional
            number of time to retry until a result is obtained
            (Default: 3)
        largest_fragment: bool, optional
            Whether to keep only the largest fragment or return full molecule
            (Default: True)

    Returns
    -------
        modified molecule or original molecule if failed
    """
    action_list = ["add", "replace"]
    new_mol = None
    new_fragment = random.choice(fragmentlist)
    new_fragment = copy.copy(new_fragment)
    action = random.choice(action_list)
    if action == "add":
        new_mol = random_fragment_add(
            mol, new_fragment, None, None, largest_fragment=largest_fragment
        )
    else:
        new_mol = random_fragment_replace(mol, new_fragment, largest_fragment=largest_fragment)

    if largest_fragment:
        return dm.keep_largest_fragment(new_mol)
    return new_mol


def random_fragment_exchange(
    mol1,
    mol2,
    largest_fragment=True,
    **kwargs,
):
    """Perform a random fragment exchange between two molecules

    Arguments
    ----------
        mol1: <Chem.Mol>
            input molecule 1
        mol2: <Chem.Mol>
            input molecule 1
        largest_fragment: bool, optional
            Whether to keep only the largest fragment or return full molecule
            (Default: True)

    Returns
    -------
        modified_mol1, modified_mol2
            Molecules obtained by exchanging fragment between mol1 and mol2.
            In case of failure, mol1, mol2 are returned

    """
    mol1_fragments = dm.fragment.anybreak(mol1, remove_parent=False)
    mol2_fragments = dm.fragment.anybreak(mol2, remove_parent=False)
    random.shuffle(mol1_fragments)
    random.shuffle(mol2_fragments)
    child_mols1 = list(
        all_fragment_attach(
            random.choice(mol1_fragments),
            mol2_fragments,
            bond_between_rings=True,
            max_num_action=10,
            asMols=True,
        )
    )
    child_mols2 = list(
        all_fragment_attach(
            random.choice(mol2_fragments),
            mol1_fragments,
            bond_between_rings=True,
            max_num_action=10,
            asMols=True,
        )
    )
    random.shuffle(child_mols1)
    random.shuffle(child_mols2)
    child_mol1 = dm.sanitize_best(child_mols1) or mol1
    child_mol2 = dm.sanitize_best(child_mols2) or mol2
    del child_mols1
    del child_mols2

    if largest_fragment:
        return dm.keep_largest_fragment(child_mol1), dm.keep_largest_fragment(child_mol2)
    return child_mol1, child_mol2


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


def random_prune_mol(
    mol,
    largest_fragment=True,
    min_size=10,
    max_size=25,
    **kwargs,
):
    """Break molecule into set of fragment using recap only (to maintain synthesis path)
    then returns one of the fragments, with probability proportional to their sizes.
    Worst case scenario, nothing happens and the original molecule will be returned.

    Arguments
    ----------
        mol: <Chem.Mol>
            Input molecule
        largest_fragment: bool, optional
            In case molecule is still fragmented, return largest fragment only
        min_size: int, optional
            minimum size of fragment, as in number of atoms
        max_size: int, optional
            maximum size of the fragment, as in number of atoms

    Returns
    -------
        Pruned molecule
    """
    breaker = kwargs.get("breaker", dm.fragment.brics)
    mass = mol.GetNumHeavyAtoms()  # Descriptors.ExactMolWt(mol)
    if mass > max_size:
        frags = breaker(mol, remove_parent=True, sanitize=True)
        frag_smiles = set([])
        for x in frags:
            if (
                x is not None
                and x.GetNumHeavyAtoms() < max_size
                and x.GetNumHeavyAtoms() > min_size
            ):
                try:
                    frag_smiles.add(dm.to_smiles(x))
                except:
                    pass

        frag_smiles = list(frag_smiles)
        frag_mols = [dm.to_mol(x) for x in frag_smiles]
        frag_weights = [len(sm) for sm in frag_smiles]
        if len(frag_mols) > 0:
            truncate = random.choices(frag_mols, k=1, weights=frag_weights)
        else:
            truncate = [None]
        mol = dm.sanitize_best(truncate) or mol
    if largest_fragment:
        return dm.keep_largest_fragment(mol)
    return mol


def random_reaction_apply(mol, rxns, rcts, mapping, **kwargs):
    """Randomly apply a chemical reaction from the reaction list, to the molecule

    Arguments
    ---------
        mol: <Chem.Mol>
            input molecule
        rxns: list
            List of reactions
        rcts: list
            List of molecules to use as reactants
        mapping: dict of dict of list
            mapping of the reactants to their position in the reactions
            e.g: {reaction_idx:{reactant_pos1:list_of_matching_blocks_idx2, reactant_pos2:list_of_matching_blocks_idx2}}

    Returns
    -------
        Molecule after random reaction application or original molecule if failed
    """

    rxn_list = [(int(k), rxns[int(k)]) for k in mapping.keys()]
    can_rxn = [dm.reactions.can_react(rxn, mol) for _, rxn in rxn_list]
    filt_rxns = [i for i, val in enumerate(can_rxn) if val >= 0]
    if filt_rxns:
        rxn_pos = random.choice(filt_rxns)
        reactants = []
        rxn_i, rxn = rxn_list[rxn_pos]
        if rxn_i not in mapping:
            rct_map = mapping[str(rxn_i)]
        else:
            rct_map = mapping[rxn_i]
        for r_pos in range(rxn.GetNumReactantTemplates()):
            tmp_ = rct_map[str(r_pos)] if r_pos not in rct_map else rct_map[r_pos]
            candidate = random.choice(tmp_)
            reactants.append(rcts[int(candidate)])
        rct_pos = can_rxn[rxn_pos]
        molist = [
            x
            for x in dm.reactions.apply_reaction(
                (rxn, tuple(reactants)), mol, rct_pos, single_output=False
            )
            if x is not None
        ]
        mol = random.choice(molist) if molist else mol
    return mol


def random_reaction_merge(mol1, mol2, rxns, **kwargs):
    """Merge two molecules by attempting to apply a reaction that matches both

    Arguments
    ----------
        mol1: <Chem.Mol>
            input molecule 1
        mol2: <Chem.Mol>
            input molecule 1
        rxns: list
            List of reactions
    Returns
    -------
        modified_mol1, modified_mol2
            Molecules obtained by exchanging fragment between mol1 and mol2.
            In case of failure, fragments are exchanged instead and in the worse
            case, mol1, mol2 are returned

    """
    mol_out_1, mol_out_2 = None, None
    rxn_list = list(enumerate(rxns))
    can_rxn = [
        (dm.reactions.can_react(rxn, mol1), dm.reactions.can_react(rxn, mol2))
        for _, rxn in rxn_list
    ]
    filt_rxns = [
        i for i, (val1, val2) in enumerate(can_rxn) if val1 >= 0 and val2 >= 0 and (val1 != val2)
    ]
    mol_outcomes = []
    while len(filt_rxns) > 0:
        random.shuffle(filt_rxns)
        rxn_pos = filt_rxns.pop(0)
        reactants = [None, None]
        rct_pos1, rct_pos2 = can_rxn[rxn_pos]
        reactants[rct_pos1] = mol1
        reactants[rct_pos2] = mol2
        _, rxn = rxn_list[rxn_pos]
        molist = rxn.RunReactants(tuple(reactants))
        mol_outcomes.extend(dm.reactions.compute_reaction_product(molist, single_output=False))

    random.shuffle(mol_outcomes)
    molist = iter(mol_outcomes)
    mol_out_1 = next(molist, None)
    mol_out_2 = next(molist, None)
    if mol_out_1 and mol_out_2:
        return mol_out_1, mol_out_2

    # do fragment exchange instead
    mol1_fragments = dm.fragment.anybreak(mol1, remove_parent=False, sanitize=True)
    mol2_fragments = dm.fragment.anybreak(mol2, remove_parent=False, sanitize=True)
    random.shuffle(mol1_fragments)
    random.shuffle(mol2_fragments)
    if len(mol1_fragments) > 0 and len(mol2_fragments) > 0:
        child_mols1 = list(
            all_fragment_update(
                random.choice(mol1_fragments),
                mol2_fragments,
                bond_between_rings=True,
                max_num_action=5,
                asMols=True,
            )
        )
        child_mols2 = list(
            all_fragment_update(
                random.choice(mol2_fragments),
                mol1_fragments,
                bond_between_rings=True,
                max_num_action=5,
                asMols=True,
            )
        )
        random.shuffle(child_mols1)
        random.shuffle(child_mols2)
        mol_out_1 = mol_out_1 or dm.sanitize_best(child_mols1)
        mol_out_2 = mol_out_2 or dm.sanitize_best(child_mols2)
        del child_mols1
        del child_mols2

    mol_out_1 = mol_out_1 or mol1
    mol_out_2 = mol_out_2 or mol2
    return dm.keep_largest_fragment(mol_out_1), dm.keep_largest_fragment(mol_out_2)
