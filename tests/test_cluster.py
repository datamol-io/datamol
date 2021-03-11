import unittest

import numpy as np

import datamol as dm


class TestCluster(unittest.TestCase):
    def test_cluster_mols(self):

        # Get some mols
        data = dm.data.freesolv()
        smiles = data["smiles"].iloc[:100].tolist()
        mols = [dm.to_mol(s) for s in smiles]

        _, mol_clusters = dm.cluster_mols(mols, cutoff=0.7)
        cluster_sizes = [15, 12, 3, 6, 9, 9, 4, 1, 4, 3, 3, 2, 3]
        assert [len(c) for c in mol_clusters[:13]] == cluster_sizes

    def test_pick_diverse(self):

        # Get some mols
        data = dm.data.freesolv()
        smiles = data["smiles"].iloc[:100].tolist()
        mols = [dm.to_mol(s) for s in smiles]

        indices, _ = dm.pick_diverse(mols, npick=18, seed=19)

        excepted_indices = np.array(
            [9, 43, 32, 89, 74, 96, 42, 91, 17, 67, 56, 94, 98, 16, 66, 58, 93, 3]
        )

        assert np.all(indices == excepted_indices)

    def test_pick_centroids(self):
        data = dm.data.freesolv()
        smiles = data["smiles"].iloc[:100].tolist()
        mols = [dm.to_mol(s) for s in smiles]
        indices, centroids = dm.pick_centroids(
            mols, npick=18, threshold=0.7, method="sphere", n_jobs=-1
        )
        excepted_indices = np.array(
            [0, 1, 2, 3, 4, 5, 8, 11, 13, 15, 16, 17, 18, 19, 21, 23, 25, 32]
        )

        assert np.all(indices == excepted_indices)

    def test_assign_to_centroids(self):
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
