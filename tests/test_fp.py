import warnings

import datamol as dm


def test_to_fp():

    smiles = "CC(=O)Oc1ccccc1C(=O)O"
    mol = dm.to_mol(smiles)

    assert dm.to_fp(mol).shape[0] == 2048
    assert dm.to_fp(mol).sum() == 29


def test_list_fp():
    assert set(dm.list_supported_fingerprints().keys()) == {
        "atompair",
        "atompair-count",
        "avalon-count",
        "ecfp",
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
    }


def test_all_fps():

    smiles = "CC(=O)Oc1ccccc1C(=O)O"
    mol = dm.to_mol(smiles)

    fp_infos = {}
    for fp_type in dm.list_supported_fingerprints():

        print(fp_type)
        args = {}
        args["mol"] = mol
        args["as_array"] = True
        args["fp_type"] = fp_type
        fp = dm.to_fp(**args)

        fp_infos[fp_type] = dict(size=len(fp), bits_sum=fp.sum())

    print(fp_infos)

    assert fp_infos == {
        "maccs": {"size": 167, "bits_sum": 21},
        "ecfp": {"size": 2048, "bits_sum": 29},
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
    }


def test_fp_deprecated_args_warnings():
    smiles = "CC(=O)Oc1ccccc1C(=O)O"
    mol = dm.to_mol(smiles)

    args = {}
    args["mol"] = mol
    args["radius"] = 3
    args["fp_size"] = 2048
    args["useFeatures"] = True
    args["as_array"] = True
    args["fp_type"] = "ecfp"

    with warnings.catch_warnings(record=True) as w:
        dm.to_fp(**args)

        assert len(w) == 1
        assert issubclass(w[-1].category, DeprecationWarning)
        assert "will be removed in datamol 0.4.0" in str(w[-1].message)

    args = {}
    args["mol"] = mol
    args["use_features"] = True
    args["as_array"] = True
    args["fp_type"] = "ecfp"

    with warnings.catch_warnings(record=True) as w:
        dm.to_fp(**args)

        assert len(w) == 1
        assert issubclass(w[-1].category, DeprecationWarning)
        assert "will be removed in datamol 0.4.0" in str(w[-1].message)
