import unittest

import datamol as dm


class TestSimilarity(unittest.TestCase):
    def test_pdist(self):

        smiles_list = ["CC(=O)Oc1ccccc1C(=O)O", "C1OC1CC", "c1cc2ccccc2cc1"]
        mols = [dm.to_mol(smiles) for smiles in smiles_list]

        dist_mat, valid_idx = dm.pdist(mols)

        assert dist_mat.shape == (3, 3)
        assert dist_mat.sum() == 5.661904761904761

        assert valid_idx.shape == (3,)
        assert valid_idx.sum() == 3

        dist_mat, valid_idx = dm.pdist(mols, n_jobs=None)

        assert dist_mat.shape == (3, 3)
        assert dist_mat.sum() == 5.661904761904761

        assert valid_idx.shape == (3,)
        assert valid_idx.sum() == 3
