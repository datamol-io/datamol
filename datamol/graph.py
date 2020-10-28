from rdkit import Chem


def _get_networkx():
    try:
        import networkx as nx

        return nx
    except ImportError:
        raise ImportError("You must install networkx from https://networkx.org/.")


def to_graph(mol: Chem.Mol):
    """Convert a molecule to a network x graph. A list of properties are added
    to every nodes and edges.

    Args:
        mol (Chem.Mol): a molecule.

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
    mol: Chem.Mol,
    atom_idx_1: int,
    atom_idx_2: int,
    ignore_cycle_basis: bool = False,
):
    """Get all simple path between two atoms of a molecule

    Args:
        mol (Chem.Mol): a molecule
        atom_idx_1 (int): Atom index 1.
        atom_idx_2 (int): Atom index 2.
        ignore_cycle_basis: Whether to ignore cycle basis.
            Defaults to False.

    Returns:
        [type]: [description]
    """

    nx = _get_networkx()

    adj = Chem.rdmolops.GetAdjacencyMatrix(mol)
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
