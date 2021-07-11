import datamol as dm


def test_to_fp():

    smiles = "CC(=O)Oc1ccccc1C(=O)O"
    mol = dm.to_mol(smiles)

    assert dm.to_fp(mol).shape[0] == 2048
    assert dm.to_fp(mol).sum() == 29


def test_list_fp():
    assert dm.list_supported_fingerprints() == {
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
        args["fp_args"] = dict()
        args["as_array"] = True
        args["fp_type"] = fp_type
        fp = dm.to_fp(**args)

        fp_infos[fp_type] = dict(size=len(fp), bits_sum=fp.sum())

    assert fp_infos == {
        "ecfp-count": {"size": 2048, "bits_sum": 42},
        "atompair": {"size": 2048, "bits_sum": 68},
        "layered": {"size": 2048, "bits_sum": 335},
        "avalon-count": {"size": 512, "bits_sum": 168},
        "fcfp-count": {"size": 2048, "bits_sum": 42},
        "pattern": {"size": 2048, "bits_sum": 173},
        "atompair-count": {"size": 2048, "bits_sum": 78},
        "maccs": {"size": 167, "bits_sum": 21},
        "rdkit": {"size": 2048, "bits_sum": 354},
        "estate": {"size": 79, "bits_sum": 13},
        "topological": {"size": 2048, "bits_sum": 18},
        "erg": {"size": 315, "bits_sum": 23.4},
        "ecfp": {"size": 2048, "bits_sum": 29},
        "topological-count": {"size": 2048, "bits_sum": 19},
    }
