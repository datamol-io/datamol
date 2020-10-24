import unittest

import datamol as dm


class TestData(unittest.TestCase):
    def test_freesolv(self):

        data = dm.data.freesolv()
        assert data.shape == (642, 4)
        assert list(data.columns) == ["iupac", "smiles", "expt", "calc"]
