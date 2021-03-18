import gc
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
    assert mol.GetConformer(0).GetPositions().shape == (4, 3)

    smiles = "CCCC"
    mol = dm.to_mol(smiles)
    mol = dm.conformers.generate(mol, rms_cutoff=None, minimize_energy=True)
    assert mol.GetNumConformers() == 50
    assert "rdkit_uff_energy" in mol.GetConformer(0).GetPropsAsDict()

    smiles = "CCCC"
    mol = dm.to_mol(smiles)
    mol = dm.conformers.generate(mol, rms_cutoff=1, minimize_energy=False)
    assert mol.GetNumConformers() == 23
    assert mol.GetConformer(0).GetPositions().shape == (4, 3)

    smiles = "CCCC"
    mol = dm.to_mol(smiles)
    mol = dm.conformers.generate(mol, rms_cutoff=1, minimize_energy=True)
    assert mol.GetNumConformers() == 25
    assert "rdkit_uff_energy" in mol.GetConformer(0).GetPropsAsDict()


@pytest.mark.skip_platform("win")
def test_sasa():

    with pytest.raises(ValueError):
        smiles = "O=C(C)Oc1ccccc1C(=O)O"
        mol = dm.to_mol(smiles)
        mol = dm.conformers.sasa(mol)

    smiles = "CCCC=O"
    mol = dm.to_mol(smiles)
    mol = dm.conformers.generate(mol, minimize_energy=True)
    sasa = dm.conformers.sasa(mol)
    assert sasa.shape == (50,)


def test_rmsd():

    with pytest.raises(ValueError):
        smiles = "O=C(C)Oc1ccccc1C(=O)O"
        mol = dm.to_mol(smiles)
        mol = dm.conformers.rmsd(mol)

    smiles = "CCCC=O"
    mol = dm.to_mol(smiles)
    mol = dm.conformers.generate(mol, rms_cutoff=None, minimize_energy=True)
    rmsd = dm.conformers.rmsd(mol)
    assert rmsd.shape == (50, 50)


def test_cluster():
    # no centroids
    smiles = "O=C(C)Oc1ccccc1C(=O)O"
    mol = dm.to_mol(smiles)
    mol = dm.conformers.generate(mol, rms_cutoff=None)
    clustered_mol = dm.conformers.cluster(mol, centroids=False)
    assert len(clustered_mol) == 2
    assert clustered_mol[0].GetNumConformers() == 40
    assert clustered_mol[1].GetNumConformers() == 10

    # centroids
    smiles = "O=C(C)Oc1ccccc1C(=O)O"
    mol = dm.to_mol(smiles)
    mol = dm.conformers.generate(mol, rms_cutoff=None)
    clustered_mol = dm.conformers.cluster(mol, centroids=True)
    assert clustered_mol.GetNumConformers() == 2

    # no centroids - minimize
    smiles = "O=C(C)Oc1ccccc1C(=O)O"
    mol = dm.to_mol(smiles)
    mol = dm.conformers.generate(mol, rms_cutoff=None, minimize_energy=True)
    clustered_mol = dm.conformers.cluster(mol, centroids=False)
    assert len(clustered_mol) == 2
    assert clustered_mol[0].GetNumConformers() == 44
    assert clustered_mol[1].GetNumConformers() == 6

    # centroids - minimize
    smiles = "O=C(C)Oc1ccccc1C(=O)O"
    mol = dm.to_mol(smiles)
    mol = dm.conformers.generate(mol, rms_cutoff=None, minimize_energy=True)
    clustered_mol = dm.conformers.cluster(mol, centroids=True)
    assert clustered_mol.GetNumConformers() == 2
