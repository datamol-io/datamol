# type: ignore
import pytest

import numpy as np
import datamol as dm
import random


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


# We disable this test as it takes too long due to the size of the SMILES
@pytest.mark.skip
def test_generate_fail():
    # This mol should fail
    smiles = "C=CC1=C(N)Oc2cc1c(-c1cc(C(C)O)cc(=O)cc1C1NCC(=O)N1)c(OC)c2OC"
    mol = dm.to_mol(smiles)
    assert dm.conformers.generate(mol, n_confs=1, verbose=True, ignore_failure=True) is None


def test_generate_2():
    smiles = "CCCC"
    mol = dm.to_mol(smiles)
    mol = dm.conformers.generate(mol, rms_cutoff=None, minimize_energy=True)
    assert mol.GetNumConformers() == 50
    assert "rdkit_UFF_energy" in mol.GetConformer(0).GetPropsAsDict()


def test_generate_3():
    smiles = "CCCC"
    mol = dm.to_mol(smiles)
    mol = dm.conformers.generate(mol, rms_cutoff=1, minimize_energy=False)
    # The exact cluster population depends on the RDKit release. The public
    # contract is that RMS pruning removes redundant conformers.
    assert 0 < mol.GetNumConformers() < 50
    assert mol.GetConformer(0).GetPositions().shape == (4, 3)


def test_generate_4():
    smiles = "CCCC"
    mol = dm.to_mol(smiles)
    mol = dm.conformers.generate(mol, rms_cutoff=1, minimize_energy=True)
    assert 0 < mol.GetNumConformers() < 50
    assert all(
        "rdkit_UFF_energy" in conformer.GetPropsAsDict() for conformer in mol.GetConformers()
    )


def test_generate_rms_cutoff_applies_to_returned_heavy_atom_conformers():
    mol = dm.to_mol("CCCCCC")
    mol = dm.conformers.generate(
        mol,
        n_confs=50,
        rms_cutoff=0.25,
        minimize_energy=True,
        align_conformers=True,
    )

    rmsd = dm.conformers.rmsd(mol)
    pairwise_rmsd = rmsd[np.triu_indices_from(rmsd, k=1)]
    assert np.all(pairwise_rmsd >= 0.25)


def test_generate_reports_unsatisfied_stereochemical_constraints():
    # This SMILES requests trans double bonds in a small ring. Modern RDKit
    # correctly rejects that impossible stereochemical geometry.
    mol = dm.to_mol("[H]/C1=C(\\[H])C([H])([H])/N=C(\\[H])OC([H])([H])C1=O")

    with pytest.raises(ValueError, match="enforce_chirality=False"):
        dm.conformers.generate(mol, n_confs=1)


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
    cluster_sizes = [cluster.GetNumConformers() for cluster in clustered_mol]
    assert sum(cluster_sizes) == 50
    assert cluster_sizes[0] > cluster_sizes[1] > 0

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
    cluster_sizes = [cluster.GetNumConformers() for cluster in clustered_mol]
    assert sum(cluster_sizes) == 50
    assert cluster_sizes[0] > cluster_sizes[1] > 0

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

    np.testing.assert_array_almost_equal(coords, new_coords - 10, decimal=0)


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
    assert len(scores) == len(mols)
    assert len(aligned_mols) == len(mols)

    for i, (mol, aligned_mol) in enumerate(zip(mols, aligned_mols)):
        p1 = mol.GetConformer().GetPositions()
        p2 = aligned_mol.GetConformer().GetPositions()

        if i == 0:
            # The first molecule is the reference so the positions should remain the same.
            assert np.allclose(p1, p2)
        else:
            assert not np.allclose(p1, p2)


def test_align_conformers_O3A():
    smiles_list = [
        "Nc1cnn(-c2ccccc2)c(=O)c1Cl",
        "Cc1ccn(-c2ccccc2)c(=O)c1F",
        "Cc1cnn(-c2ccccc2)c(=O)c1Cl",
        "Cc1cnn(-c2ccccc2)c(=O)c1",
    ]
    mols = [dm.to_mol(smiles) for smiles in smiles_list]
    mols = [dm.conformers.generate(mol, n_confs=1) for mol in mols]

    # Align
    aligned_mols, scores = dm.conformers.align_conformers(mols, backend="O3A")

    # Check
    assert len(scores) == len(mols)
    assert len(aligned_mols) == len(mols)

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


def test_align_conformers_wrong_backend():
    smiles_list = [
        "Nc1cnn(-c2ccccc2)c(=O)c1Cl",
        "Cc1ccn(-c2ccccc2)c(=O)c1F",
        "Cc1cnn(-c2ccccc2)c(=O)c1Cl",
        "Cc1cnn(-c2ccccc2)c(=O)c1",
    ]
    mols = [dm.to_mol(smiles) for smiles in smiles_list]
    mols = [dm.conformers.generate(mol, n_confs=1) for mol in mols]

    with pytest.raises(ValueError):
        dm.conformers.align_conformers(mols, backend="not_supported")


def test_conformers_minimized_sorted():
    smiles = "CC(=O)NC1CCC2=CC(=C(C(=C2C3=CC=C(C(=O)C=C13)OC)OC)OC)OC"
    mol = dm.to_mol(smiles)
    mol = dm.conformers.generate(mol, n_confs=10, minimize_energy=True)

    confs = list(mol.GetConformers())
    energies = np.array([conf.GetPropsAsDict()["rdkit_UFF_energy"] for conf in confs])

    assert np.all(np.diff(energies) >= 0)


def test_conformers_non_minimized_sorted():
    smiles = "CC(=O)NC1CCC2=CC(=C(C(=C2C3=CC=C(C(=O)C=C13)OC)OC)OC)OC"
    mol = dm.to_mol(smiles)
    mol = dm.conformers.generate(mol, n_confs=10, minimize_energy=False)

    confs = list(mol.GetConformers())
    energies = np.array([conf.GetPropsAsDict()["rdkit_UFF_energy"] for conf in confs])

    assert np.all(np.diff(energies) >= 0)


def test_keep_conformers_from_indice():
    mol = dm.to_mol("CC")
    mol = dm.conformers.generate(mol, n_confs=10)

    mol2 = dm.conformers.keep_conformers(mol, 5)
    conf_ids = [conf.GetId() for conf in mol2.GetConformers()]

    assert mol2.GetNumConformers() == 1
    assert conf_ids == [0]


def test_keep_conformers_from_indice_default_conf():
    mol = dm.to_mol("CC")
    mol = dm.conformers.generate(mol, n_confs=10)

    mol2 = dm.conformers.keep_conformers(mol)
    conf_ids = [conf.GetId() for conf in mol2.GetConformers()]

    assert mol2.GetNumConformers() == 1
    assert conf_ids == [0]


def test_keep_conformers_from_indices():
    mol = dm.to_mol("CC")
    mol = dm.conformers.generate(mol, n_confs=10)

    mol2 = dm.conformers.keep_conformers(mol, [0, 4, 6])
    conf_ids = [conf.GetId() for conf in mol2.GetConformers()]

    assert mol2.GetNumConformers() == 3
    assert conf_ids == [0, 1, 2]


def test_keep_conformers_from_indices_keep_ids():
    mol = dm.to_mol("CC")
    mol = dm.conformers.generate(mol, n_confs=10)

    mol2 = dm.conformers.keep_conformers(mol, [0, 4, 6], assign_id=False)
    conf_ids = [conf.GetId() for conf in mol2.GetConformers()]

    assert mol2.GetNumConformers() == 3
    assert conf_ids == [0, 4, 6]


def test_conformer_energy():
    mol = dm.to_mol("O=C(C)Oc1ccccc1C(=O)O")

    random.seed(42)
    np.random.seed(42)
    mol1 = dm.conformers.generate(mol, ewindow=7)
    assert 0 < mol1.GetNumConformers() <= 50
    uff_props = [conf.GetPropsAsDict() for conf in mol1.GetConformers()]
    uff_energies = [props["rdkit_UFF_energy"] for props in uff_props]
    uff_deltas = [props["rdkit_UFF_delta_energy"] for props in uff_props]
    assert uff_energies == sorted(uff_energies)
    assert uff_deltas[0] == pytest.approx(0)
    assert max(uff_deltas) <= 7

    mol2 = dm.conformers.generate(mol, forcefield="MMFF94s", eratio=3)
    assert 0 < mol2.GetNumConformers() <= 50
    mmff_props = [conf.GetPropsAsDict() for conf in mol2.GetConformers()]
    mmff_energies = [props["rdkit_MMFF94s_energy"] for props in mmff_props]
    mmff_deltas = [props["rdkit_MMFF94s_delta_energy"] for props in mmff_props]
    assert mmff_energies == sorted(mmff_energies)
    assert mmff_deltas[0] == pytest.approx(0)
    assert max(mmff_deltas) <= 3 * dm.descriptors.n_rotatable_bonds(dm.add_hs(mol))

    mol3 = dm.conformers.generate(mol, forcefield="MMFF94s_noEstat", minimize_energy=True)
    no_estat_props = [conf.GetPropsAsDict() for conf in mol3.GetConformers()]
    no_estat_energies = [props["rdkit_MMFF94s_noEstat_energy"] for props in no_estat_props]
    no_estat_deltas = [props["rdkit_MMFF94s_noEstat_delta_energy"] for props in no_estat_props]
    assert no_estat_energies == sorted(no_estat_energies)
    assert no_estat_deltas[0] == pytest.approx(0)


def test_conformer_no_rotatable_bonds():
    mol = dm.to_mol("c1ccccc1")

    random.seed(42)
    np.random.seed(42)
    mol1 = dm.conformers.generate(mol, minimize_energy=True)

    random.seed(42)
    np.random.seed(42)
    mol2 = dm.conformers.generate(mol, minimize_energy=True, eratio=3)

    # `eratio` should be ignored for this molecule as it has no rotatable bonds
    assert mol1.GetNumConformers() == mol2.GetNumConformers()
