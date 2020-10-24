import unittest
import pathlib

import datamol as dm

ROOT_DIR = pathlib.Path(__file__).parent.resolve()
DATA_DIR = ROOT_DIR / "data"


class TestMol(unittest.TestCase):
    def test_to_mol(self):
        smiles = "O=C(C)Oc1ccccc1C(=O)O"
        mol = dm.to_mol(smiles)
        assert mol.GetNumAtoms() == 13
