import pytest

import datamol as dm


def test_to_fp():
    smiles = "CC(=O)Oc1ccccc1C(=O)O"
    mol = dm.to_mol(smiles)

    assert dm.to_fp(mol).shape[0] == 2048
    assert dm.to_fp(mol).sum() == 31


def test_list_fp():
    assert set(dm.list_supported_fingerprints().keys()) == {
        "atompair",
        "atompair-count",
        "avalon-count",
        "ecfp",
        "fcfp",
        "ecfp-count",
        "erg",
        "estate",
        "fcfp-count",
        "layered",
        "maccs",
        "pattern",
        "rdkit",
        "topological",
        "topological-count",
        "rdkit-count",
    }


def test_all_fps():
    smiles = "CC(=O)Oc1ccccc1C(=O)O"
    mol = dm.to_mol(smiles)

    fp_infos = {}
    for fp_type in dm.list_supported_fingerprints():
        fold_size = None
        if fp_type == "rdkit-count":
            fold_size = 2048

        print(fp_type)
        args = {}
        args["mol"] = mol
        args["as_array"] = True
        args["fp_type"] = fp_type
        args["fold_size"] = fold_size
        fp = dm.to_fp(**args)

        fp_infos[fp_type] = dict(size=len(fp), bits_sum=fp.sum())

    print(fp_infos)

    assert fp_infos == {
        "maccs": {"size": 167, "bits_sum": 21},
        "ecfp": {"size": 2048, "bits_sum": 31},
        "fcfp": {"size": 2048, "bits_sum": 22},
        "topological": {"size": 2048, "bits_sum": 18},
        "atompair": {"size": 2048, "bits_sum": 68},
        "rdkit": {"size": 2048, "bits_sum": 354},
        "pattern": {"size": 2048, "bits_sum": 173},
        "layered": {"size": 2048, "bits_sum": 335},
        "erg": {"size": 315, "bits_sum": 23.4},
        "estate": {"size": 79, "bits_sum": 13},
        "avalon-count": {"size": 512, "bits_sum": 168},
        "ecfp-count": {"size": 2048, "bits_sum": 35},
        "fcfp-count": {"size": 2048, "bits_sum": 35},
        "topological-count": {"size": 2048, "bits_sum": 19},
        "atompair-count": {"size": 2048, "bits_sum": 78},
        "rdkit-count": {"size": 2048, "bits_sum": 301},
    }


def test_fp_invalid_input():
    args = {}
    args["mol"] = None
    args["radius"] = 3

    with pytest.raises(ValueError):
        dm.to_fp(**args)

    args["mol"] = "dsdsdsd"
    with pytest.raises(ValueError):
        dm.to_fp(**args)
