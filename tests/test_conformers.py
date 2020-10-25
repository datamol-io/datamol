import unittest
import pytest

import datamol as dm


class TestConformers(unittest.TestCase):
    def test_generate(self):

        with pytest.raises(ValueError):
            smiles = "CCCC"
            mol = dm.to_mol(smiles)
            mol = dm.conformers.generate(mol, method="custom_method")

        smiles = "CCCC"
        mol = dm.to_mol(smiles)
        mol = dm.conformers.generate(mol)
        assert mol.GetNumConformers() == 50

        conf = mol.GetConformer(0)
        assert conf.GetPositions().shape == (14, 3)

        props = conf.GetPropsAsDict()
        assert "rdkit_uff_energy" in props

    def test_sasa(self):

        with pytest.raises(ValueError):
            smiles = "O=C(C)Oc1ccccc1C(=O)O"
            mol = dm.to_mol(smiles)
            mol = dm.conformers.sasa(mol)

        smiles = "CCCC=O"
        mol = dm.to_mol(smiles)
        mol = dm.conformers.generate(mol)
        sasa = dm.conformers.sasa(mol)
        assert sasa.shape == (50,)

    def test_rmsd(self):

        with pytest.raises(ValueError):
            smiles = "O=C(C)Oc1ccccc1C(=O)O"
            mol = dm.to_mol(smiles)
            mol = dm.conformers.rmsd(mol)

        smiles = "CCCC=O"
        mol = dm.to_mol(smiles)
        mol = dm.conformers.generate(mol)
        rmsd = dm.conformers.rmsd(mol)
        assert rmsd.shape == (50, 50)
