import pytest

import pandas as pd
from rdkit import Chem

import datamol as dm


def test_to_smiles():

    smiles = "O=C(C)Oc1ccccc1C(=O)O"
    mol = dm.to_mol(smiles)

    smiles = dm.to_smiles(
        mol,
        isomeric=True,
        ordered=True,
        explicit_bonds=False,
        explicit_hs=False,
    )
    assert smiles == "CC(=O)Oc1ccccc1C(=O)O"

    smiles = dm.to_smiles(
        mol,
        isomeric=True,
        ordered=False,
        explicit_bonds=True,
        explicit_hs=False,
    )
    assert smiles == "C-C(=O)-O-c1:c:c:c:c:c:1-C(=O)-O"

    smiles = dm.to_smiles(
        mol,
        isomeric=True,
        ordered=False,
        explicit_bonds=False,
        explicit_hs=True,
    )
    assert smiles == "[CH3][C](=[O])[O][c]1[cH][cH][cH][cH][c]1[C](=[O])[OH]"

    smiles = "O=C(C)Oc1ccccc1C(=O)O"
    mol = dm.to_mol(smiles)
    randomized_smiles = dm.to_smiles(mol, randomize=True)
    randomized_mol = dm.to_mol(randomized_smiles)

    assert dm.to_smiles(randomized_mol) == dm.to_smiles(mol)


def test_to_selfies():
    smiles = "CC(=O)Oc1ccccc1C(=O)O"
    mol = dm.to_mol(smiles)

    true_sf = (
        "[C][C][Branch1_2][C][=O][O][C][=C][C][=C][C][=C][Ring1][Branch1_2][C][Branch1_2][C][=O][O]"
    )

    selfies = dm.to_selfies(smiles)
    assert selfies == true_sf

    selfies = dm.to_selfies(mol)
    assert selfies == true_sf


def test_from_selfies():
    selfies = (
        "[C][C][Branch1_2][C][=O][O][C][=C][C][=C][C][=C][Ring1][Branch1_2][C][Branch1_2][C][=O][O]"
    )

    smiles = dm.from_selfies(selfies, as_mol=False)
    assert smiles == "CC(=O)OC1=CC=CC=C1C(=O)O"

    mol = dm.from_selfies(selfies, as_mol=True)
    assert dm.to_smiles(mol) == "CC(=O)Oc1ccccc1C(=O)O"


def test_to_smarts():
    smiles = "O=C(C)Oc1ccccc1C(=O)O"
    mol = dm.to_mol(smiles)

    smarts = dm.to_smarts(mol, keep_hs=True)
    assert smarts == "[CH3]-[C](=[O])-[O]-[c]1:[cH]:[cH]:[cH]:[cH]:[c]:1-[C](=[O])-[OH]"

    smarts = dm.to_smarts(mol, keep_hs=False)
    assert smarts == "[CH3]-[C](=[O])-[O]-[c]1:[cH]:[cH]:[cH]:[cH]:[c]:1-[C](=[O])-[OH]"

    assert dm.to_smarts(None) is None


def test_inchi():
    smiles = "CC(=O)Oc1ccccc1C(=O)O"
    mol = dm.to_mol(smiles)

    inchi = dm.to_inchi(mol)
    assert inchi == "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)"

    inchikey = dm.to_inchikey(mol)
    assert inchikey == "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"

    new_mol = dm.from_inchi(inchi)
    assert dm.to_smiles(new_mol) == smiles

    assert dm.to_inchi(None) is None
    assert dm.to_inchikey(None) is None
    assert dm.from_inchi(None) is None


def test_to_df(datadir):

    data_path = datadir / "TUBB3-observations.sdf"
    mols = dm.read_sdf(data_path)
    df = dm.to_df(mols)

    assert df.shape == (10, 12)
    assert list(df.columns) == [
        "smiles",
        "zinc_id",
        "ortholog_name",
        "gene_name",
        "affinity",
        "chembldocid",
        "title",
        "reference.pubmed_id",
        "reference.doi",
        "reference.chembl_id",
        "reference.journal",
        "reference.year",
    ]


def test_from_df(datadir):
    data_path = datadir / "TUBB3-observations.sdf"
    df = dm.read_sdf(data_path, as_df=True)

    mols = dm.from_df(df)

    assert len(mols) == 10
    assert isinstance(mols[0], Chem.Mol)

    assert set(mols[0].GetPropsAsDict().keys()) == {
        "zinc_id",
        "ortholog_name",
        "gene_name",
        "affinity",
        "chembldocid",
        "title",
        "reference.pubmed_id",
        "reference.doi",
        "reference.chembl_id",
        "reference.journal",
        "reference.year",
    }

    assert dm.from_df(pd.DataFrame()) == []


def test_from_df_conserve_smiles(datadir):
    data_path = datadir / "freesolv.csv"
    df = dm.read_csv(data_path)
    mols = dm.from_df(df, conserve_smiles=True)
    assert "smiles" in mols[0].GetPropsAsDict().keys()


def test_to_df_smiles_warning(datadir, caplog):
    data_path = datadir / "freesolv.csv"
    df = dm.read_csv(data_path)

    mols = dm.from_df(df, conserve_smiles=True)
    df = dm.to_df(mols)

    assert sum(df.columns == "smiles") == 2

    for record in caplog.records:
        assert record.levelname != "WARNING"
    assert (
        "The SMILES column name provided ('smiles') is already present in the properties of the molecules"
        not in caplog.text
    )


def test_to_cxsmiles():
    mol = dm.to_mol("OC1=CC2CCCCC2[N:1]=C1")
    smiles = dm.to_smiles(mol, cxsmiles=True)
    assert smiles == "OC1=CC2CCCCC2[N:1]=C1 |atomProp:9.molAtomMapNumber.1|"


def test_to_smiles_fail():
    smiles = dm.to_smiles(55, allow_to_fail=False)
    assert smiles == None

    # NOTE(hadim): ideally you want to catch only `Boost.Python.ArgumentError` here.
    with pytest.raises(Exception):
        dm.to_smiles(55, allow_to_fail=True)
