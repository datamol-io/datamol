from rdkit.Chem.rdmolops import GetAdjacencyMatrix
from rdkit.Chem.AllChem import RenumberAtoms
from typing import Dict, List, Union
from loguru import logger
from datamol.mol import add_hs
from datamol import Mol


def _get_networkx():
    try:
        import networkx as nx

        return nx
    except ImportError:
        raise ImportError("You must install networkx from https://networkx.org/.")


def to_graph(mol: Mol):
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
    mol: Mol,
    atom_idx_1: int,
    atom_idx_2: int,
    ignore_cycle_basis: bool = False,
):
    """Get all simple path between two atoms of a molecule

    Args:
        mol (dm.Mol): a molecule
        atom_idx_1 (int): Atom index 1.
        atom_idx_2 (int): Atom index 2.
        ignore_cycle_basis: Whether to ignore cycle basis.
            Defaults to False.

    Returns:
        [type]: [description]
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
    mol1: Mol,
    mol2: Mol,
    match_atoms_on: List[str] = ["atomic_num"],
    match_bonds_on: List[str] = ["bond_type"],
    explicit_hs: bool = False,
) -> List[Dict[int, int]]:
    """
    Match the node indices of 2 molecular graphs, with optional usage of atomic number and edge type.

    Args:
        mol1, mol2: The molecules to match their indices.
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

        explicit_hs: Whether to consider the hydrogens explicitly when matching.

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

    # Add explicit hydrogens
    if explicit_hs:
        mol1 = add_hs(mol1)
        mol2 = add_hs(mol2)

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
    mol: Mol,
    mol_template: Mol,
    enforce_atomic_num: bool = False,
    enforce_bond_type: bool = False,
    explicit_hs: bool = False,
    verbose: bool = True,
) -> Union[Mol, type(None)]:
    """
    Re-order the nodes of a molecular graph from the nodes of a template molecule.
    Molecular graphs and atom types need to be identical, but edge types and charges
    are not enforced.

    This is particularily useful when dealing with XYZ files containing node ordering,
    but with missing information regarding charges and edge types.

    If you only need to match bond orders, you can check the function
    `rdkit.Chem.AllChem.AssignBondOrdersFromTemplate`.

    Args:
        mol: The molecule to re-order
        mol_template: The molecule containing the right node order.
        enforce_atomic_num: Whether to enforce atomic number. Atomic numbers are always enforced
            for a first try. If no match are found and this parameter is `False`,
            the matching is tried again.
        enforce_bond_type: Whether to enforce bond types. Bond types are always enforced
            for a first try. If no match are found and this parameter is `False`,
            the matching is tried again.
        explicit_hs: Whether to consider the hydrogens explicitly when matching.
        verbose: Whether to warn when the matching does not work

    Returns:
        - `None` if the molecular graphs do not match (both the graph and atom types).
            Pring a warning.
        - `None` if multiple matche are found, which can happen for symmetric molecules such as benzene
            Pring a warning.
        - `Mol` The re-ordered molecule when a single match is found.
    """

    # Match the ordering of the graphs
    matches = match_molecular_graphs(
        mol_template,
        mol,
        match_atoms_on=["atomic_num"],
        match_bonds_on=["bond_type"],
        explicit_hs=explicit_hs,
    )

    # If no matches were found, retry without bond types
    if (len(matches) == 0) and (not enforce_bond_type):
        matches = match_molecular_graphs(
            mol_template,
            mol,
            match_atoms_on=["atomic_num"],
            match_bonds_on=[],
            explicit_hs=explicit_hs,
        )

    # If no matches were found, retry without atom types
    if (len(matches) == 0) and (not enforce_atomic_num):
        matches = match_molecular_graphs(
            mol_template,
            mol,
            match_atoms_on=[],
            match_bonds_on=["bond_type"],
            explicit_hs=explicit_hs,
        )

    # If no matches were found, retry without bond and atom types
    if (len(matches) == 0) and (not enforce_bond_type) and (not enforce_atomic_num):
        matches = match_molecular_graphs(
            mol_template, mol, match_atoms_on=[], match_bonds_on=[], explicit_hs=explicit_hs
        )

    # If no match were found, exit the function and return None
    if len(matches) == 0:
        if verbose:
            logger.warning("No match was found")
        return None

    # If many matches were found, exit the function and return None
    if len(matches) > 1:
        if verbose:
            logger.warning(f"{len(matches)} matches were found, ordering is ambiguous")
        return None

    # Re-order the molecule from the matching indices of the template
    match = matches[0]
    match = [match[ii] for ii in range(mol.GetNumAtoms())]
    reordered_mol = RenumberAtoms(mol, match)

    return reordered_mol
