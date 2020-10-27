import unittest

import numpy as np

import datamol as dm


class TestFragment(unittest.TestCase):
    def test_brics(self):
        smiles = "CCCOCc1cc(c2ncccc2)ccc1"
        mol = dm.to_mol(smiles)
        frags = dm.fragment.brics(mol)
        assert len(frags) == 9

    def test_frag(self):
        smiles = "CCCOCc1cc(c2ncccc2)ccc1"
        mol = dm.to_mol(smiles)
        frags = dm.fragment.frag(mol)
        assert len(frags) == 9

    def test_recap(self):
        smiles = "CCCOCc1cc(c2ncccc2)ccc1"
        mol = dm.to_mol(smiles)
        frags = dm.fragment.recap(mol)
        assert len(frags) == 3

    def test_anybreak(self):
        smiles = "CCCOCc1cc(c2ncccc2)ccc1"
        mol = dm.to_mol(smiles)
        frags = dm.fragment.anybreak(mol)
        assert len(frags) == 9

    def test_mmpa(self):
        smiles = "CCCOCc1cc(c2ncccc2)ccc1"
        mol = dm.to_mol(smiles)

        frags = dm.fragment.mmpa_cut(mol)
        assert len(frags) == 39
        assert (
            "CCCOCc1cccc(-c2ccccn2)c1,C(C[*:2])[*:1],C[*:1].c1ccc(-c2cccc(CO[*:2])c2)nc1\n" in frags
        )

    def test_assemble(self):

        # Fragment a molecule
        smiles = "CCCOCc1cc(c2ncccc2)ccc1"
        mol = dm.to_mol(smiles)
        frags = dm.fragment.brics(mol)

        # Limit the number of fragments to work with because
        # assembling is computationally intensive.
        frags = frags[:2]

        # Assemble molecules from the list of fragments
        mols = list(dm.fragment.assemble_fragment_order(frags, max_n_mols=4))

        assert len(mols) == 4
