import datamol as dm


def test_to_mol():

    smiles = "O=C(C)Oc1ccccc1C(=O)O"
    mol = dm.to_mol(smiles)
    graph = dm.to_graph(mol)

    assert list(graph.nodes) == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    assert list(graph.edges) == [
        (0, 1),
        (1, 2),
        (1, 3),
        (3, 4),
        (4, 5),
        (4, 9),
        (5, 6),
        (6, 7),
        (7, 8),
        (8, 9),
        (9, 10),
        (10, 11),
        (10, 12),
    ]


def test_get_all_path_between():
    smiles = "c1cc2cccccc2c1"
    mol = dm.to_mol(smiles)

    all_paths = dm.get_all_path_between(mol, 8, 4, ignore_cycle_basis=False)
    assert all_paths == [[8, 2, 3, 4], [8, 7, 6, 5, 4], [8, 9, 0, 1, 2, 3, 4]]

    all_paths = dm.get_all_path_between(mol, 8, 4, ignore_cycle_basis=True)
    assert all_paths == [[8, 2, 3, 4], [8, 7, 6, 5, 4]]
