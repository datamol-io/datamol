# type: ignore
import pytest

import numpy as np
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

    mol = dm.to_mol(smiles)
    mol = dm.conformers.generate(
        mol, rms_cutoff=None, minimize_energy=False, embed_params={"useRandomCoords": True}
    )
    assert mol.GetNumConformers() == 50
    assert mol.GetConformer(0).GetPositions().shape == (4, 3)

    # This mol should fail
    smiles = "C=CC1=C(N)Oc2cc1c(-c1cc(C(C)O)cc(=O)cc1C1NCC(=O)N1)c(OC)c2OC"
    mol = dm.to_mol(smiles)
    assert dm.conformers.generate(mol, n_confs=1, verbose=True, ignore_failure=True) is None

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
    assert clustered_mol[0].GetNumConformers() > 30
    assert clustered_mol[1].GetNumConformers() > 5

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
    assert clustered_mol[0].GetNumConformers() > 30
    assert clustered_mol[1].GetNumConformers() > 5

    # centroids - minimize
    smiles = "O=C(C)Oc1ccccc1C(=O)O"
    mol = dm.to_mol(smiles)
    mol = dm.conformers.generate(mol, rms_cutoff=None, minimize_energy=True)
    clustered_mol = dm.conformers.cluster(mol, centroids=True)
    assert clustered_mol.GetNumConformers() == 2


def test_get_coords():
    mol = dm.to_mol("CC")
    mol = dm.conformers.generate(mol, n_confs=1)

    assert dm.conformers.get_coords(mol).shape == (2, 3)


def test_center_of_mass():
    mol = dm.to_mol("CC")
    mol = dm.conformers.generate(mol, n_confs=1)

    # geomtric center
    center = dm.conformers.center_of_mass(mol, use_atoms=False)
    coords = dm.conformers.get_coords(mol)
    np.testing.assert_array_almost_equal(coords.mean(axis=0), center)

    # mass center
    center = dm.conformers.center_of_mass(mol, use_atoms=True)
    assert center.shape == (3,)


def test_translate():
    mol = dm.to_mol("CC")
    mol = dm.conformers.generate(mol, n_confs=1)

    coords = dm.conformers.get_coords(mol)
    print(coords)

    dm.conformers.translate(mol, [10, 10, 10])
    new_coords = dm.conformers.get_coords(mol)
    print(new_coords - 10)

    np.testing.assert_array_almost_equal(coords, new_coords - 10, decimal=1)


def test_align_conformers():
    smiles_list = [
        "Nc1cnn(-c2ccccc2)c(=O)c1Cl",
        "Cc1ccn(-c2ccccc2)c(=O)c1F",
        "Cc1cnn(-c2ccccc2)c(=O)c1Cl",
        "Cc1cnn(-c2ccccc2)c(=O)c1",
    ]
    mols = [dm.to_mol(smiles) for smiles in smiles_list]
    mols = [dm.conformers.generate(mol, n_confs=1) for mol in mols]

    # Align
    aligned_mols, scores = dm.conformers.align_conformers(mols)

    # Check
    assert scores.shape == (4,)

    for i, (mol, aligned_mol) in enumerate(zip(mols, aligned_mols)):
        p1 = mol.GetConformer().GetPositions()
        p2 = aligned_mol.GetConformer().GetPositions()

        if i == 0:
            # The first molecule is the reference so the positions should remain the same.
            assert np.allclose(p1, p2)
        else:
            assert not np.allclose(p1, p2)


def test_align_conformers_without_conformer():
    smiles_list = [
        "Nc1cnn(-c2ccccc2)c(=O)c1Cl",
        "Cc1ccn(-c2ccccc2)c(=O)c1F",
        "Cc1cnn(-c2ccccc2)c(=O)c1Cl",
        "Cc1cnn(-c2ccccc2)c(=O)c1",
    ]
    mols = [dm.to_mol(smiles) for smiles in smiles_list]

    with pytest.raises(ValueError):
        dm.conformers.align_conformers(mols)
