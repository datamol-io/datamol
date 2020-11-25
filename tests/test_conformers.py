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
        mol = dm.conformers.generate(mol, n_confs=2, minimize_energy=True)
        assert mol.GetNumConformers() == 2

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
        mol = dm.conformers.generate(mol, minimize_energy=False)
        sasa = dm.conformers.sasa(mol)
        assert sasa.shape == (50,)

    def test_rmsd(self):

        with pytest.raises(ValueError):
            smiles = "O=C(C)Oc1ccccc1C(=O)O"
            mol = dm.to_mol(smiles)
            mol = dm.conformers.rmsd(mol)

        smiles = "CCCC=O"
        mol = dm.to_mol(smiles)
        mol = dm.conformers.generate(mol, minimize_energy=False)
        rmsd = dm.conformers.rmsd(mol)
        assert rmsd.shape == (50, 50)

    def test_cluster(self):

        smiles = "O=C(C)Oc1ccccc1C(=O)O"
        mol = dm.to_mol(smiles)
        mol = dm.conformers.generate(mol)
        mol.GetNumConformers()

        clustered_mol = dm.conformers.cluster(mol, return_centroids=True)
        assert clustered_mol.GetNumConformers() == 3
