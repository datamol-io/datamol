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
