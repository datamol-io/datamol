import pytest

import datamol as dm


def test_generate():

    with pytest.raises(ValueError):
        smiles = "CCCC"
        mol = dm.to_mol(smiles)
        mol = dm.conformers.generate(mol, method="custom_method")

    smiles = "CCCC"
    mol = dm.to_mol(smiles)
    mol = dm.conformers.generate(mol, rms_cutoff=None, minimize_energy=False)
    assert mol.GetNumConformers() == 50

    conf = mol.GetConformer(0)
    assert conf.GetPositions().shape == (14, 3)

    smiles = "CCCC"
    mol = dm.to_mol(smiles)
    mol = dm.conformers.generate(mol, rms_cutoff=1, minimize_energy=False)
    assert mol.GetNumConformers() == 23

    # NOTE(hadim): `minimize_energy=True` fails on GA.
    # props = conf.GetPropsAsDict()
    # assert "rdkit_uff_energy" in props


def test_sasa():

    with pytest.raises(ValueError):
        smiles = "O=C(C)Oc1ccccc1C(=O)O"
        mol = dm.to_mol(smiles)
        mol = dm.conformers.sasa(mol)

    smiles = "CCCC=O"
    mol = dm.to_mol(smiles)
    mol = dm.conformers.generate(mol, minimize_energy=False)
    sasa = dm.conformers.sasa(mol)
    assert sasa.shape == (17,)


def test_rmsd():

    with pytest.raises(ValueError):
        smiles = "O=C(C)Oc1ccccc1C(=O)O"
        mol = dm.to_mol(smiles)
        mol = dm.conformers.rmsd(mol)

    smiles = "CCCC=O"
    mol = dm.to_mol(smiles)
    mol = dm.conformers.generate(mol, minimize_energy=False)
    rmsd = dm.conformers.rmsd(mol)
    assert rmsd.shape == (17, 17)


def test_cluster():
    # no centroids
    smiles = "O=C(C)Oc1ccccc1C(=O)O"
    mol = dm.to_mol(smiles)
    mol = dm.conformers.generate(mol, n_confs=10, rms_cutoff=None, minimize_energy=False)
    clustered_mol = dm.conformers.cluster(mol, centroids=False)
    assert len(clustered_mol) == 3
    assert clustered_mol[0].GetNumConformers() == 5
    assert clustered_mol[1].GetNumConformers() == 3
    assert clustered_mol[2].GetNumConformers() == 2

    # centroids
    smiles = "O=C(C)Oc1ccccc1C(=O)O"
    mol = dm.to_mol(smiles)
    mol = dm.conformers.generate(mol, n_confs=10, rms_cutoff=None, minimize_energy=False)
    clustered_mol = dm.conformers.cluster(mol, centroids=True)
    assert clustered_mol.GetNumConformers() == 3
