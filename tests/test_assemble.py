import unittest

import numpy as np

import datamol as dm


class TestAssemble(unittest.TestCase):
    def test_assemble_brics(self):

        # Fragment a molecule
        smiles = "CCCOCc1cc(c2ncccc2)ccc1"
        mol = dm.to_mol(smiles)
        frags = dm.fragment.brics(mol)

        # Limit the number of fragments to work with because
        # assembling is computationally intensive.
        frags = frags[:2]

        # Assemble molecules from the list of fragments
        mols = list(dm.assemble.assemble_brics_order(frags, max_n_mols=4))

        assert len(mols) == 4
