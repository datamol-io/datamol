import unittest
import pytest

import datamol as dm

N_CONFS = 2


def test_generate():

    with pytest.raises(ValueError):
        smiles = "CCCC"
        mol = dm.to_mol(smiles)
        mol = dm.conformers.generate(mol, n_confs=10, method="custom_method")

    smiles = "CCCC"
    mol = dm.to_mol(smiles)
    mol = dm.conformers.generate(mol, n_confs=N_CONFS)
    assert mol.GetNumConformers() == N_CONFS

    conf = mol.GetConformer(0)
    assert conf.GetPositions().shape == (14, 3)

    props = conf.GetPropsAsDict()
    assert "rdkit_uff_energy" in props


def test_sasa():

    with pytest.raises(ValueError):
        smiles = "O=C(C)Oc1ccccc1C(=O)O"
        mol = dm.to_mol(smiles)
        mol = dm.conformers.sasa(mol)

    smiles = "CCCC=O"
    mol = dm.to_mol(smiles)
    mol = dm.conformers.generate(mol, n_confs=N_CONFS)
    sasa = dm.conformers.sasa(mol)
    assert sasa.shape == (N_CONFS,)


def test_rmsd():

    with pytest.raises(ValueError):
        smiles = "O=C(C)Oc1ccccc1C(=O)O"
        mol = dm.to_mol(smiles)
        mol = dm.conformers.rmsd(mol)

    smiles = "CCCC=O"
    mol = dm.to_mol(smiles)
    mol = dm.conformers.generate(mol, n_confs=N_CONFS)
    rmsd = dm.conformers.rmsd(mol)
    assert rmsd.shape == (N_CONFS, N_CONFS)


def test_cluster():

    smiles = "O=C(C)Oc1ccccc1C(=O)O"
    mol = dm.to_mol(smiles)
    mol = dm.conformers.generate(mol, n_confs=N_CONFS)
    mol.GetNumConformers()

    clustered_mol = dm.conformers.cluster(mol, return_centroids=True)
    assert clustered_mol.GetNumConformers() == 2
