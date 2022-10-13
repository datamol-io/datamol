from typing import Dict
from typing import List
from typing import Union
from typing import Optional

from rdkit.Chem.rdmolops import GetAdjacencyMatrix
from rdkit.Chem.AllChem import RenumberAtoms  # type: ignore
from loguru import logger

import datamol as dm


def _get_networkx():
    try:
        import networkx as nx

        return nx
    except ImportError:
        raise ImportError("You must install networkx from https://networkx.org/.")


def to_graph(mol: dm.Mol):
    """Convert a molecule to a network x graph. A list of properties are added
    to every nodes and edges.

    Args:
        mol (dm.Mol): a molecule.

    Returns:
        mol_graph (networkx.Graph): a graph representing the molecule.
    """

    nx = _get_networkx()

    mol_graph = nx.Graph()
    for atom in mol.GetAtoms():
        mol_graph.add_node(
            atom.GetIdx(),
            atomic_num=atom.GetAtomicNum(),
            formal_charge=atom.GetFormalCharge(),
            chiral_tag=atom.GetChiralTag(),
            hybridization=atom.GetHybridization(),
            num_explicit_hs=atom.GetNumExplicitHs(),
            implicit_valence=atom.GetImplicitValence(),
            degree=atom.GetDegree(),
            symbol=atom.GetSymbol(),
            ring_atom=atom.IsInRing(),
            is_aromatic=atom.GetIsAromatic(),
        )
    for bond in mol.GetBonds():
        mol_graph.add_edge(
            bond.GetBeginAtomIdx(),
            bond.GetEndAtomIdx(),
            bond_type=bond.GetBondType(),
        )
    return mol_graph


def get_all_path_between(
    mol: dm.Mol,
    atom_idx_1: int,
    atom_idx_2: int,
    ignore_cycle_basis: bool = False,
) -> list:
    """Get all simple path between two atoms of a molecule

    Args:
        mol (dm.Mol): a molecule
        atom_idx_1 (int): Atom index 1.
        atom_idx_2 (int): Atom index 2.
        ignore_cycle_basis: Whether to ignore cycle basis.
            Defaults to False.

    Returns:
        list of path between two atoms.
    """

    nx = _get_networkx()

    adj = GetAdjacencyMatrix(mol)
    G = nx.Graph(adj)
    path = nx.all_simple_paths(G, source=atom_idx_1, target=atom_idx_2)

    if ignore_cycle_basis:
        rings = [set(x) for x in mol.GetRingInfo().AtomRings()]
        final_path = []
        for p in path:
            reject_path = False
            for r in rings:
                if r.issubset(set(p)):
                    reject_path = True
                    break
            if not reject_path:
                final_path.append(p)
        path = final_path

    return list(path)


def match_molecular_graphs(
    mol1: dm.Mol,
    mol2: dm.Mol,
    match_atoms_on: List[str] = ["atomic_num"],
    match_bonds_on: List[str] = ["bond_type"],
) -> List[Dict[int, int]]:
    """
    Match the node indices of 2 molecular graphs,
    with optional usage of atomic number and edge type.

    Note:
        The matching fails if the hydrogens are implicit in one molecule,
        but explicit in the other.

    Note:
        Explicit hydrogens might lead to too many matches, since for an atom with 2
        hydrogens, they can be re-ordered in any way.

    Args:
        mol1: A molecule.
        mol2: A molecule.
        match_atoms_on: Properties on which to match the atom types.
            By default, it matches on the `'atomic_num'` property.
            Empty list means that it does not consider atom features during matching.

            Other properties are defined by the `datamol.graph.to_graph` function and include:
            - atomic_num
            - formal_charge
            - chiral_tag
            - hybridization
            - num_explicit_hs
            - implicit_valence
            - degree
            - symbol
            - ring_atom
            - is_aromatic

        match_bonds_on: Properties on which to match the bond types.
            Empty list means that it does not consider bond features during matching.
            By default, it matches on the `'bond_type'` property.
            No other properties are defined by the `datamol.graph.to_graph` function.

    Returns:
        A list of all matches dictionaries. In case of a single match, the list has len==1.
        Each dictionary contains as key the indices of `mol1` and as value the corresponding
        indices of `mol2`.
    """

    nx = _get_networkx()

    if isinstance(match_atoms_on, str):
        match_atoms_on = [match_atoms_on]
    if isinstance(match_bonds_on, str):
        match_bonds_on = [match_bonds_on]

    def node_match_fn(node1, node2):
        """Function that matches the atomic number"""
        return all([node1[prop] == node2[prop] for prop in match_atoms_on])

    def edge_match_fn(edge1, edge2):
        """Function that matches the bond type"""
        return all([edge1[prop] == edge2[prop] for prop in match_bonds_on])

    # Convert to networkx graph
    g1 = to_graph(mol1)
    g2 = to_graph(mol2)

    # Use the `match` function to find the matching indices
    node_match = node_match_fn if len(match_atoms_on) > 0 else None
    edge_match = edge_match_fn if len(match_bonds_on) > 0 else None
    graph_matcher = nx.algorithms.isomorphism.vf2userfunc.GraphMatcher(
        g1, g2, node_match=node_match, edge_match=edge_match
    )
    matches = list(graph_matcher.match())

    return matches


def reorder_mol_from_template(
    mol: dm.Mol,
    mol_template: dm.Mol,
    enforce_atomic_num: bool = False,
    enforce_bond_type: bool = False,
    ambiguous_match_mode: str = "No",
    verbose: bool = True,
) -> Optional[dm.Mol]:
    """
    Re-order the nodes of a molecular graph from the nodes of a template molecule.
    Molecular graphs and atom types need to be identical, but edge types and charges
    are not enforced.

    This is particularily useful when dealing with XYZ files containing node ordering,
    but with missing information regarding charges and edge types.

    !!! note

        * If you only need to match bond orders, you can check the function
        `rdkit.Chem.AllChem.AssignBondOrdersFromTemplate`.
        * The matching fails if the hydrogens are implicit in one molecule,
        but explicit in the other.
        * Explicit hydrogens might lead to too many matches, since for an atom with 2
        hydrogens, they can be re-ordered in any way.

    Args:
        mol: The molecule to re-order
        mol_template: The molecule containing the right node order.
        enforce_atomic_num: Whether to enforce atomic number. Atomic numbers are always enforced
            for a first try. If no match are found and this parameter is `False`,
            the matching is tried again.
        enforce_bond_type: Whether to enforce bond types. Bond types are always enforced
            for a first try. If no match are found and this parameter is `False`,
            the matching is tried again.
        ambiguous_match_mode: Whether to allow ambiguous matching. This means that,
            if there are many matches to the molecule, it will still re-order
            the molecule according to specific rules. Options are:
            - "no": Does not allow ambiguous matching.
            - "hs-only": Allow matching of ambiguous hydrogens. Does not work if trying
              to match implicit with explicit hydrogens.
            - "first": Return the first match.
            - "best": Return the match with the least errors on atom type, edges type, and edge stereo.
              Errors on the atoms are counted with 1 point, on the charge with 0.25 points,
              on the edges with 0.25 points, and on the Stereo with 0.05 points.
              If the option `enforce_atomic_num` is used, then no errors on the atoms are allowed.
              If the option `enforce_bond_type` is used, then no errors on the edges are allowed.
            - "best-first": "best", followed by "first".
        verbose: Whether to warn when the matching does not work or is ambiguous.
            Different warnings are raised depending on the value of `ambiguous_match_mode`.

    Returns:
        - `None` if the molecular graphs do not match (both the graph and atom types).
            Pring a warning.
        - `None` if multiple matche are found, which can happen for symmetric molecules such as benzene
            Pring a warning.
        - `Mol` The re-ordered molecule when a single match is found.
    """

    ambiguous_match_mode = ambiguous_match_mode.lower()

    # Match the ordering of the graphs
    matches = match_molecular_graphs(
        mol_template,
        mol,
        match_atoms_on=["atomic_num"],
        match_bonds_on=["bond_type"],
    )

    # If no matches were found, retry without bond types
    if (len(matches) == 0) and (not enforce_bond_type):
        matches = match_molecular_graphs(
            mol_template,
            mol,
            match_atoms_on=["atomic_num"],
            match_bonds_on=[],
        )

    # If no matches were found, retry without atom types
    if (len(matches) == 0) and (not enforce_atomic_num):
        matches = match_molecular_graphs(
            mol_template,
            mol,
            match_atoms_on=[],
            match_bonds_on=["bond_type"],
        )

    # If no matches were found, retry without bond and atom types
    if (len(matches) == 0) and (not enforce_bond_type) and (not enforce_atomic_num):
        matches = match_molecular_graphs(mol_template, mol, match_atoms_on=[], match_bonds_on=[])

    # If no match were found, exit the function and return None
    if len(matches) == 0:
        if verbose:
            logger.warning("No match was found")
        return None

    if len(matches) > 1:
        # In case we want to allow ambiguous match of hydrogens
        if ambiguous_match_mode == "hs-only":
            first_keys = list(matches[0].keys())
            all_hs_mismatch = True
            for this_match in matches:
                this_keys = list(this_match.keys())
                keys_mismatch = [
                    ii for ii in range(len(this_keys)) if (first_keys[ii] != this_keys[ii])
                ]
                atoms_mismatch = [mol.GetAtomWithIdx(key).GetAtomicNum() for key in keys_mismatch]
                all_hs = all([atom == 1 for atom in atoms_mismatch])
                if not all_hs:
                    all_hs_mismatch = False
                    break
            if all_hs_mismatch:
                matches = matches[0:1]
            else:
                if verbose:
                    logger.warning(
                        f"{len(matches)} matches were found, ordering is ambiguous, even when ignoring hydrogens"
                    )
                return None

        # Compute the number of atoms and bonds mismatch, and select the one with the least mismatch
        if (ambiguous_match_mode in ["best", "best-first"]) and not (
            enforce_atomic_num and enforce_bond_type
        ):
            num_mismatches = []
            for this_match in matches:
                num_atoms_mismatch, num_charge_mismatch = 0, 0

                # Get the number of atomic mismatch
                for key, val in this_match.items():
                    atom1 = mol.GetAtomWithIdx(val)
                    atom2 = mol_template.GetAtomWithIdx(key)
                    num_atoms_mismatch += atom1.GetAtomicNum() != atom2.GetAtomicNum()
                    num_charge_mismatch += atom1.GetFormalCharge() != atom2.GetFormalCharge()

                # Get the number of bond mismatch
                num_bonds_type_mismatch, num_bonds_stereo_mismatch = 0, 0
                for bond1 in mol_template.GetBonds():
                    begin_idx, end_idx = bond1.GetBeginAtomIdx(), bond1.GetEndAtomIdx()
                    bond2 = mol.GetBondBetweenAtoms(this_match[begin_idx], this_match[end_idx])
                    num_bonds_type_mismatch += bond1.GetBondType() != bond2.GetBondType()
                    num_bonds_stereo_mismatch += (bond1.GetStereo() != bond2.GetStereo()) or (
                        bond1.GetBondDir() != bond2.GetBondDir()
                    )

                num_mismatches.append(
                    (1 * num_atoms_mismatch)
                    + (0.25 * num_charge_mismatch)
                    + (0.25 * num_bonds_type_mismatch)
                    + (0.05 * num_bonds_stereo_mismatch)
                )
            min_mismatch_idx = [
                ii for ii in range(len(num_mismatches)) if num_mismatches[ii] == min(num_mismatches)
            ]
            matches = [matches[idx] for idx in min_mismatch_idx]

        # Select the first matching element
        if ambiguous_match_mode in ["first", "best-first"]:
            matches = [matches[0]]

    if len(matches) > 1:
        # If many matches were found, exit the function and return None
        if ambiguous_match_mode == "no":
            if verbose:
                logger.warning(f"{len(matches)} matches were found, ordering is ambiguous")
            return None

    # Re-order the molecule from the matching indices of the template
    match = matches[0]
    match = [match[ii] for ii in range(mol.GetNumAtoms())]
    reordered_mol = RenumberAtoms(mol, match)

    return reordered_mol
