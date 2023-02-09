import pytest

import numpy as np
import datamol as dm
import datamol.utils.testing


def test_pdist():
    smiles_list = ["CC(=O)Oc1ccccc1C(=O)O", "C1OC1CC", "c1cc2ccccc2cc1"]
    mols = [dm.to_mol(smiles) for smiles in smiles_list]

    dist_mat = dm.pdist(mols)

    assert dist_mat.shape == (3, 3)
    assert dist_mat.sum() == 5.6757105943152455

    dist_mat = dm.pdist(mols, n_jobs=None)

    assert dist_mat.shape == (3, 3)
    assert dist_mat.sum() == 5.6757105943152455


def test_pdist_condensed():
    smiles_list = ["CC(=O)Oc1ccccc1C(=O)O", "C1OC1CC", "c1cc2ccccc2cc1"]
    mols = [dm.to_mol(smiles) for smiles in smiles_list]

    dist_mat = dm.pdist(mols, squareform=False)

    assert dist_mat.shape == (3,)
    assert dist_mat.sum() == 2.8378552971576227


def test_cdist():
    smiles_list1 = ["CC(=O)Oc1ccccc1C(=O)O", "C1OC1CC", "c1cc2ccccc2cc1"]
    mols1 = [dm.to_mol(smiles) for smiles in smiles_list1]

    smiles_list2 = [
        "COc1cc(Nc2ncc(Cl)c(-c3cccc(CC#N)c3)n2)ccc1N1CCN(C)CC1",
        "ON=C(O)CCCCCN=C(O)C=C1c2ccccc2-c2ccccc21",
        "COc1ccc(CCc2nnc(-c3ccc4nc[nH]c4c3)o2)cc1Cl",
    ]
    mols2 = [dm.to_mol(smiles) for smiles in smiles_list2]

    dist_mat = dm.cdist(mols1, mols2)

    assert dist_mat.shape == (3, 3)
    assert np.isclose(dist_mat.mean(), 0.9416270180919872)


def test_cdist_chunked():
    smiles_list1 = ["CC(=O)Oc1ccccc1C(=O)O", "C1OC1CC", "c1cc2ccccc2cc1"]
    mols1 = [dm.to_mol(smiles) for smiles in smiles_list1]

    smiles_list2 = [
        "COc1cc(Nc2ncc(Cl)c(-c3cccc(CC#N)c3)n2)ccc1N1CCN(C)CC1",
        "ON=C(O)CCCCCN=C(O)C=C1c2ccccc2-c2ccccc21",
        "COc1ccc(CCc2nnc(-c3ccc4nc[nH]c4c3)o2)cc1Cl",
    ]
    mols2 = [dm.to_mol(smiles) for smiles in smiles_list2]

    d1 = dm.cdist(mols1, mols2, distances_chunk=True)
    d2 = dm.cdist(mols1, mols2, distances_chunk=False)

    assert d1.shape == d2.shape
    assert np.allclose(d1, d2)


def test_cdist_pdist_consistent():
    smiles_list1 = ["CC(=O)Oc1ccccc1C(=O)O", "C1OC1CC", "c1cc2ccccc2cc1"]
    mols1 = [dm.to_mol(smiles) for smiles in smiles_list1]

    dist_mat = dm.cdist(mols1, mols1)
    dist_mat2 = dm.pdist(mols1)

    assert np.isclose(dist_mat.mean(), dist_mat2.mean())
    assert np.allclose(dist_mat, dist_mat2)


def test_cdist_pdist_invalid_input():
    smiles_list = ["CC(=O)Oc1ccccc1C(=O)O", "C1OC1CC", "c1cc2ccccc2cc1", "dsdsdsd"]

    with pytest.raises(ValueError):
        dm.similarity.cdist(smiles_list, smiles_list)

    with pytest.raises(ValueError):
        dm.similarity.pdist(smiles_list)


def test_datamol_pdist_same_as_rdkit():
    smiles_list = [
        "COc1cc(Nc2ncc(Cl)c(-c3cccc(CC#N)c3)n2)ccc1N1CCN(C)CC1",
        "ON=C(O)CCCCCN=C(O)C=C1c2ccccc2-c2ccccc21",
        "COc1ccc(CCc2nnc(-c3ccc4nc[nH]c4c3)o2)cc1Cl",
    ]

    dist_mat = dm.similarity.pdist(smiles_list)
    dist_mat_rdkit = datamol.utils.testing.pdist_rdkit(smiles_list)

    assert np.allclose(dist_mat, dist_mat_rdkit)


def test_datamol_cdist_same_as_rdkit():
    smiles_list = [
        "COc1cc(Nc2ncc(Cl)c(-c3cccc(CC#N)c3)n2)ccc1N1CCN(C)CC1",
        "ON=C(O)CCCCCN=C(O)C=C1c2ccccc2-c2ccccc21",
        "COc1ccc(CCc2nnc(-c3ccc4nc[nH]c4c3)o2)cc1Cl",
    ]

    smiles_list2 = ["CC(=O)Oc1ccccc1C(=O)O", "C1OC1CC", "c1cc2ccccc2cc1"]

    dist_mat = dm.similarity.cdist(smiles_list, smiles_list2)
    dist_mat_rdkit = datamol.utils.testing.cdist_rdkit(smiles_list, smiles_list2)

    assert np.allclose(dist_mat, dist_mat_rdkit)
