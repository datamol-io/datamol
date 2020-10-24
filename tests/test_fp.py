import unittest

import datamol as dm


class TestFP(unittest.TestCase):
    def test_to_fp(self):

        smiles = "CC(=O)Oc1ccccc1C(=O)O"
        mol = dm.to_mol(smiles)

        assert dm.to_fp(mol).shape[0] == 2048
        assert dm.to_fp(mol).sum() == 29
