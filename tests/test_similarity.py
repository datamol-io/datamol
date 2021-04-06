import unittest

import numpy as np

import datamol as dm


def test_pdist():

    smiles_list = ["CC(=O)Oc1ccccc1C(=O)O", "C1OC1CC", "c1cc2ccccc2cc1"]
    mols = [dm.to_mol(smiles) for smiles in smiles_list]

    dist_mat, valid_idx = dm.pdist(mols)

    assert dist_mat.shape == (3, 3)
    assert dist_mat.sum() == 5.661904761904761

    assert valid_idx.shape == (3,)
    assert valid_idx.sum() == 3

    dist_mat, valid_idx = dm.pdist(mols, n_jobs=None)

    assert dist_mat.shape == (3, 3)
    assert dist_mat.sum() == 5.661904761904761

    assert valid_idx.shape == (3,)
    assert valid_idx.sum() == 3


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
    assert np.isclose(dist_mat.mean(), 0.9416646866643121)

    # Check results are the same with `pdist`
    dist_mat = dm.cdist(mols1, mols1)
    dist_mat2, _ = dm.pdist(mols1)
    assert np.isclose(dist_mat.mean(), dist_mat2.mean())
