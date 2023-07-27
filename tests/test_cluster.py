import numpy as np

import datamol as dm


def test_cluster_mols():
    # Get some mols
    data = dm.data.freesolv()
    smiles = data["smiles"].iloc[:100].tolist()
    mols = [dm.to_mol(s) for s in smiles]

    _, mol_clusters = dm.cluster_mols(mols, cutoff=0.7)
    cluster_sizes = [11, 7, 5, 3, 3, 3, 2, 3, 2, 1, 2, 2, 1]
    assert [len(c) for c in mol_clusters[:13]] == cluster_sizes


def test_pick_diverse():
    # Get some mols
    data = dm.data.freesolv()
    smiles = data["smiles"].iloc[:100].tolist()
    mols = [dm.to_mol(s) for s in smiles]

    indices, _ = dm.pick_diverse(mols, npick=18, seed=19)

    excepted_indices = np.array(
        [9, 14, 47, 50, 56, 61, 67, 89, 83, 90, 94, 10, 0, 96, 15, 58, 71, 21]
    )

    assert np.all(indices == excepted_indices)


def test_pick_centroids():
    data = dm.data.freesolv()
    smiles = data["smiles"].iloc[:100].tolist()
    mols = [dm.to_mol(s) for s in smiles]
    indices, centroids = dm.pick_centroids(
        mols, npick=18, threshold=0.7, method="sphere", n_jobs=-1
    )
    excepted_indices = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 13, 14, 15, 16, 17, 18, 20])

    assert np.all(indices == excepted_indices)


def test_assign_to_centroids():
    data = dm.data.freesolv()
    smiles = data["smiles"].iloc[:100].tolist()
    mols = [dm.to_mol(s) for s in smiles]
    indices, centroids = dm.pick_centroids(
        mols, npick=18, threshold=0.7, method="sphere", n_jobs=-1
    )

    cluster_map, cluster_list = dm.assign_to_centroids(mols, centroids, n_jobs=-1)
    # expect centroid to be in centroid list
    assert indices[0] in cluster_map[0]
    # expect no intersection after assignment
    map_intersection = set.intersection(*map(set, cluster_map.values()))
    assert len(map_intersection) == 0
    # expect some similar molecule in a given cluster
    # assert 33 in cluster_map[0]
