import unittest
import pytest

import datamol as dm


@pytest.mark.serial
def test_generate(self):

    with pytest.raises(ValueError):
        smiles = "CCCC"
        mol = dm.to_mol(smiles)
        mol = dm.conformers.generate(mol, method="custom_method")

    smiles = "CCCC"
    mol = dm.to_mol(smiles)
    mol = dm.conformers.generate(mol, minimize_energy=False)
    assert mol.GetNumConformers() == 50

    conf = mol.GetConformer(0)
    assert conf.GetPositions().shape == (14, 3)

    # NOTE(hadim): `minimize_energy=True` fails on GA.
    # props = conf.GetPropsAsDict()
    # assert "rdkit_uff_energy" in props


@pytest.mark.serial
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


@pytest.mark.serial
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


@pytest.mark.serial
def test_cluster(self):
    # NOTE(hadim): disable here since something is wrong when running on CI.
    smiles = "O=C(C)Oc1ccccc1C(=O)O"
    mol = dm.to_mol(smiles)
    mol = dm.conformers.generate(mol, n_confs=5, minimize_energy=False)

    clustered_mol = dm.conformers.cluster(mol, return_centroids=False)
    assert len(clustered_mol) == 2
