from rdkit.Chem import AllChem
from typing import Dict, List, Union
import networkx as nx
from networkx.algorithms.isomorphism.vf2userfunc import GraphMatcher
from loguru import logger


def mol_to_nx_graph(mol: AllChem.rdchem.Mol) -> nx.Graph:
    """
    Convert a molecule to a networkx graph.

    Parameters:
        mol: Molecule to convert
    Returns:
        A networkx graph containing an `atomic_num` attribute on the nodes
        and a `weight` attribute on the edges, corresponding to
        single=1, double=2, triple=3, quadruple=4, quintuple=5, aromatic=1.5.
    """
    atomic_numx_property = {
        ii: {"atomic_num": atom.GetAtomicNum()} for ii, atom in enumerate(mol.GetAtoms())
    }
    adj = AllChem.GetAdjacencyMatrix(mol, useBO=True)
    g = nx.from_numpy_matrix(adj)
    nx.set_node_attributes(g, atomic_numx_property)
    return g


def match_molecular_graphs(
    mol1: AllChem.rdchem.Mol,
    mol2: AllChem.rdchem.Mol,
    match_atomic_num: bool = True,
    match_bond_type: bool = True,
) -> List[Dict[int, int]]:
    """
    Match the node indices of 2 molecular graphs, with optional usage of atomic number and edge type.

    Parameters:
        mol1, mol2: The molecules to match their indices.
        match_atomic_num: Whether the atomic number should be used when matching the graphs.
        match_bond_type: Whether the bond type should be used when matching the graphs.
            Bond-type is limited to single, double, triple, quadruple, quituple, and aromatic bonds.

    Returns:
        A list of all matches dictionaries. In case of a single match, the list has len==1.
        Each dictionary contains as key the indices of `mol1` and as value the corresponding
        indices of `mol2`.
    """

    def node_match_fn(node1, node2):
        """Function that matches the atomic number"""
        return node1["atomic_num"] == node2["atomic_num"]

    def edge_match_fn(edge1, edge2):
        """Function that matches the edge weights"""
        return edge1["weight"] == edge2["weight"]

    # Convert to networkx graph
    g1 = mol_to_nx_graph(mol1)
    g2 = mol_to_nx_graph(mol2)

    # Use the `match` function to find the matching indices
    node_match = node_match_fn if match_atomic_num else None
    edge_match = edge_match_fn if match_bond_type else None
    graph_matcher = GraphMatcher(g1, g2, node_match=node_match, edge_match=edge_match)
    matches = list(graph_matcher.match())

    return matches


def reorder_mol_from_template(
    mol: AllChem.rdchem.Mol,
    mol_template: AllChem.rdchem.Mol,
    enforce_atomic_num: bool = False,
    enforce_bond_type: bool = False,
) -> Union[AllChem.rdchem.Mol, type(None)]:
    """
    Re-order the nodes of a molecular graph from the nodes of a template molecule.
    Molecular graphs and atom types need to be identical, but edge types and charges
    are not enforced.

    This is particularily useful when dealing with XYZ files containing node ordering,
    but with missing information regarding charges and edge types.

    Parameters:
        mol: The molecule to re-order
        mol_template: The molecule containing the right node order.
        enforce_atomic_num: Whether to enforce atomic number. Atomic numbers are always enforced
            for a first try. If no match are found and this parameter is `False`,
            the matching is tried again.
        enforce_bond_type: Whether to enforce bond types. Bond types are always enforced
            for a first try. If no match are found and this parameter is `False`,
            the matching is tried again.

    Returns:
        - `None` if the molecular graphs do not match (both the graph and atom types).
            Pring a warning.
        - `None` if multiple matche are found, which can happen for symmetric molecules such as benzene
            Pring a warning.
        - `Mol` The re-ordered molecule when a single match is found.
    """

    # Match the ordering of the graphs
    matches = match_molecular_graphs(mol_template, mol, match_atomic_num=True, match_edge_type=True)

    # If no matches were found, retry without bond types
    if (len(matches) == 0) and (not enforce_bond_type):
        matches = match_molecular_graphs(
            mol_template, mol, match_atomic_num=True, match_edge_type=False
        )

    # If no matches were found, retry without atom types
    if (len(matches) == 0) and (not enforce_atomic_num):
        matches = match_molecular_graphs(
            mol_template, mol, match_atomic_num=False, match_edge_type=True
        )

    # If no matches were found, retry without bond and atom types
    if (len(matches) == 0) and (not enforce_bond_type) and (not enforce_atomic_num):
        matches = match_molecular_graphs(
            mol_template, mol, match_atomic_num=False, match_edge_type=False
        )

    # If no match were found, exit the function and return None
    if len(matches) == 0:
        logger.warning("No match was found")
        return None

    # If many matches were found, exit the function and return None
    if len(matches) > 1:
        logger.warning(f"{len(matches)} matches were found, ordering is ambiguous")
        print(matches)
        return None

    # Re-order the molecule from the matching indices of the template
    match = matches[0]
    match = [match[ii] for ii in range(mol.GetNumAtoms())]
    reordered_mol = AllChem.RenumberAtoms(mol, match)

    return reordered_mol
