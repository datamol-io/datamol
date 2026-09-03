import pytest

import datamol as dm
from rdkit import DataStructs


def test_to_fp():
    smiles = "CC(=O)Oc1ccccc1C(=O)O"
    mol = dm.to_mol(smiles)

    assert dm.to_fp(mol).shape[0] == 2048
    assert dm.to_fp(mol).sum() == 31


def test_list_fp():
    assert set(dm.list_supported_fingerprints().keys()) == {
        "atompair",
        "atompair-count",
        "avalon",
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

    expected = {
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
        "avalon": {"size": 512, "bits_sum": 54},
        "avalon-count": {"size": 512, "bits_sum": 168},
        "ecfp-count": {"size": 2048, "bits_sum": 42},
        "fcfp-count": {"size": 2048, "bits_sum": 35},
        "topological-count": {"size": 2048, "bits_sum": 19},
        "atompair-count": {"size": 2048, "bits_sum": 78},
        "rdkit-count": {"size": 2048, "bits_sum": 301},
    }
    assert fp_infos.keys() == expected.keys()
    for fp_type, info in fp_infos.items():
        assert info["size"] == expected[fp_type]["size"]
        assert info["bits_sum"] == pytest.approx(expected[fp_type]["bits_sum"])


def test_fp_invalid_input():
    args = {}
    args["mol"] = None
    args["radius"] = 3

    with pytest.raises(ValueError):
        dm.to_fp(**args)

    args["mol"] = "dsdsdsd"
    with pytest.raises(ValueError):
        dm.to_fp(**args)


def test_fold_count_fp_supports_explicit_bit_vectors():
    fingerprint = DataStructs.ExplicitBitVect(16)
    fingerprint.SetBitsFromList([1, 5, 9])

    folded = dm.fold_count_fp(fingerprint, dim=4)
    binary = dm.fold_count_fp(fingerprint, dim=4, binary=True)

    assert folded.tolist() == [0, 3, 0, 0]
    assert binary.tolist() == [0, 1, 0, 0]


@pytest.mark.parametrize(
    "fingerprint",
    [DataStructs.ExplicitBitVect(16), DataStructs.UIntSparseIntVect(16)],
)
def test_fold_count_fp_supports_empty_fingerprints(fingerprint):
    folded = dm.fold_count_fp(fingerprint, dim=4)

    assert folded.dtype.kind == "i"
    assert folded.tolist() == [0, 0, 0, 0]
