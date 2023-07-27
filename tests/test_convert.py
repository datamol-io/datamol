from typing import cast
import pytest

import numpy as np
import pandas as pd
from rdkit import Chem
from selfies import __version__ as selfies_version
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

    mol = dm.to_mol(mol, kekulize=True)
    smiles = dm.to_smiles(
        mol, isomeric=True, ordered=False, explicit_bonds=False, explicit_hs=False, kekulize=True
    )
    assert smiles == "CC(=O)OC1=CC=CC=C1C(=O)O"

    smiles = "O=C(C)Oc1ccccc1C(=O)O"
    mol = dm.to_mol(smiles)
    randomized_smiles = dm.to_smiles(mol, randomize=True)
    randomized_mol = dm.to_mol(randomized_smiles)

    assert dm.to_smiles(randomized_mol) == dm.to_smiles(mol)


def test_to_selfies():
    smiles = "CC(=O)Oc1ccccc1C(=O)O"
    mol = dm.to_mol(smiles)

    true_sf = (
        "[C][C][=Branch1][C][=O][O][C][=C][C][=C][C][=C][Ring1][=Branch1][C][=Branch1][C][=O][O]"
    )
    if selfies_version < "2.0.0":
        true_sf = "[C][C][Branch1_2][C][=O][O][C][=C][C][=C][C][=C][Ring1][Branch1_2][C][Branch1_2][C][=O][O]"

    selfies = dm.to_selfies(smiles)
    assert selfies == true_sf

    selfies = dm.to_selfies(mol)
    assert selfies == true_sf


def test_from_selfies():
    # work with both selfies version
    selfies = (
        "[C][C][=Branch1][C][=O][O][C][=C][C][=C][C][=C][Ring1][=Branch1][C][=Branch1][C][=O][O]"
    )
    if selfies_version < "2.0.0":
        selfies = "[C][C][Branch1_2][C][=O][O][C][=C][C][=C][C][=C][Ring1][Branch1_2][C][Branch1_2][C][=O][O]"

    smiles = dm.from_selfies(selfies, as_mol=False)
    assert smiles == "CC(=O)OC1=CC=CC=C1C(=O)O"

    mol = dm.from_selfies(selfies, as_mol=True)
    assert dm.to_smiles(mol) == "CC(=O)Oc1ccccc1C(=O)O"


def test_smiles_as_smarts():
    smiles = "O=C(C)Oc1ccccc1C(=O)O"
    mol = dm.to_mol(smiles)

    smarts = dm.smiles_as_smarts(mol, keep_hs=True)
    assert smarts == "[CH3]-[C](=[O])-[O]-[c]1:[cH]:[cH]:[cH]:[cH]:[c]:1-[C](=[O])-[OH]"

    smarts = dm.smiles_as_smarts(smiles, keep_hs=False)
    assert smarts == "[CH3]-[C](=[O])-[O]-[c]1:[cH]:[cH]:[cH]:[cH]:[c]:1-[C](=[O])-[OH]"

    assert dm.smiles_as_smarts(None) is None


def test_to_smarts():
    smarts = "[OX2H][CX3]=[OX1]"
    mol = Chem.MolFromSmarts(smarts)
    expected_out = dm.to_smarts(mol)
    assert expected_out == "[O&X2&H1][C&X3]=[O&X1]"


def test_from_smarts():
    smarts = "[OX2H][CX3]=[OX1]"
    smiles = "O=C(C)Oc1ccccc1C(=O)O"
    query = dm.from_smarts(smarts)
    mol = dm.to_mol(smiles)
    assert mol.HasSubstructMatch(query)


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
    assert dm.to_inchi("") is None
    assert dm.to_inchi("C(C)(C)(C)(C)(C)") is None
    assert dm.to_inchikey(None) is None
    assert dm.to_inchikey("") is None
    assert dm.to_inchikey("C(C)(C)(C)(C)(C)") is None
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


def test_to_df_parallel(datadir):
    data_path = datadir / "TUBB3-observations.sdf"
    mols = dm.read_sdf(data_path)

    # check parallel or sequential `to_df()` gives the same output df
    large_mol_set = np.random.choice(mols, 50)
    df_sequential = dm.to_df(large_mol_set, n_jobs=1)
    df_parallel = dm.to_df(large_mol_set, n_jobs=-1)

    pd.testing.assert_frame_equal(df_sequential, df_parallel)


def test_from_df(datadir):
    data_path = datadir / "TUBB3-observations.sdf"
    df = dm.read_sdf(data_path, as_df=True)

    mols = dm.from_df(df)

    assert len(mols) == 10
    assert isinstance(mols[0], Chem.rdchem.Mol)

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

    assert "WARNING" in caplog.text
    assert (
        "The SMILES column name provided ('smiles') is already present in the properties of the molecules"
        in caplog.text
    )


def test_to_cxsmiles():
    mol = dm.to_mol("OC1=CC2CCCCC2[N:1]=C1")
    smiles = dm.to_smiles(mol, cxsmiles=True)
    assert smiles == "OC1=CC2CCCCC2[N:1]=C1 |atomProp:9.molAtomMapNumber.1|"


def test_to_smiles_fail():
    smiles = dm.to_smiles(55, allow_to_fail=False)
    assert smiles is None

    # NOTE(hadim): ideally you want to catch only `Boost.Python.ArgumentError` here.
    with pytest.raises(Exception):
        dm.to_smiles(55, allow_to_fail=True)


def test_from_df_pop_mol_column():
    df = dm.data.freesolv().iloc[:10]  # type: ignore
    mols = [dm.to_mol(smiles) for smiles in df["smiles"]]

    df: pd.DataFrame = dm.to_df(mols, mol_column="mol")  # type: ignore
    df["dummy"] = "hello"

    # test with provided mol column
    mols = dm.from_df(df.copy(), mol_column="mol")
    assert set(mols[0].GetPropsAsDict().keys()) == {"smiles", "dummy"}

    # test with automatic mol column detection
    mols = dm.from_df(df.copy())
    assert set(mols[0].GetPropsAsDict().keys()) == {"smiles", "dummy"}


def test_non_standard_inchi():
    mol1 = dm.to_mol("N=C(N)O")
    mol2 = dm.to_mol("NC(N)=O")
    mol3 = dm.to_mol("c1ccccc1")
    mol4 = dm.to_mol("c1cccc([*])c1")
    # with a standard inchi, both molecules should be identical
    assert dm.to_inchi(mol1) == "InChI=1S/CH4N2O/c2-1(3)4/h(H4,2,3,4)"
    assert dm.to_inchi(mol2) == "InChI=1S/CH4N2O/c2-1(3)4/h(H4,2,3,4)"
    assert dm.to_inchikey(mol1) == "XSQUKJJJFZCRTK-UHFFFAOYSA-N"
    assert dm.to_inchikey(mol2) == "XSQUKJJJFZCRTK-UHFFFAOYSA-N"

    # with additional layers, inchi values should be different
    assert dm.to_inchi_non_standard(mol1) == "InChI=1/CH4N2O/c2-1(3)4/h(H4,2,3,4)/f/h2,4H,3H2/b2-1?"
    assert dm.to_inchi_non_standard(mol2) == "InChI=1/CH4N2O/c2-1(3)4/h(H4,2,3,4)/f/h2-3H2"
    assert dm.to_inchikey_non_standard(mol1) == "XSQUKJJJFZCRTK-ZIALIONUNA-N"
    assert dm.to_inchikey_non_standard(mol2) == "XSQUKJJJFZCRTK-UBUOBULFNA-N"

    # for molecule without movile hydrogens the inchi and inchikey will still be different
    # because of the absence of the letter `S` in the prefix
    assert dm.to_inchi_non_standard(mol3) == "InChI=1/C6H6/c1-2-4-6-5-3-1/h1-6H"
    assert dm.to_inchi(mol3) == "InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H"
    assert dm.to_inchikey_non_standard(mol3) == "UHOVQNZJYSORNB-UHFFFAOYNA-N"
    assert dm.to_inchikey(mol3) == "UHOVQNZJYSORNB-UHFFFAOYSA-N"

    # let's check we can reload the correct molecule from the non standard inchi
    # The SMILES is unlikely to change with a simple molecule like that
    mol1_reloaded = dm.from_inchi(dm.to_inchi_non_standard(mol1))
    mol2_reloaded = dm.from_inchi(dm.to_inchi_non_standard(mol2))

    assert dm.to_smiles(mol1_reloaded) == "N=C(N)O"
    assert dm.to_smiles(mol2_reloaded) == "NC(N)=O"

    # inchikey computation should return None on query
    assert dm.to_inchikey(mol4) is None


def test_non_standard_inchi_with_options():
    mol1 = dm.to_mol("C[C@@H]1CCCC[C@@H]1O")
    mol2 = dm.to_mol("C[C@H]1CCCC[C@@H]1O")

    assert dm.to_inchikey(mol1) == "NDVWOBYBJYUSMF-RQJHMYQMSA-N"
    assert dm.to_inchikey_non_standard(mol1) == "NDVWOBYBJYUSMF-RQJHMYQMNA-N"
    assert dm.to_inchikey_non_standard(mol1, options=["/SRel"]) == "NDVWOBYBJYUSMF-JHPDDGAFNA-N"

    assert dm.to_inchikey(mol2) == "NDVWOBYBJYUSMF-BQBZGAKWSA-N"
    assert dm.to_inchikey_non_standard(mol2) == "NDVWOBYBJYUSMF-BQBZGAKWNA-N"
    assert dm.to_inchikey_non_standard(mol2, options=["/SRel"]) == "NDVWOBYBJYUSMF-WZTWBHKBNA-N"


def test_non_standard_inchi_edge_cases():
    assert dm.to_inchi_non_standard(None) is None
    assert dm.to_inchikey_non_standard(None) is None

    # test with wrong query molecules
    assert dm.to_inchi_non_standard("NC([*])C(N)=O") is None
    assert dm.to_inchikey_non_standard("NC([*])C(N)=O") is None

    # test with empty molecule
    assert dm.to_inchi_non_standard("CCC(C)(C)(C)(C)(C)") is None
    assert dm.to_inchikey_non_standard("CCC(C)(C)(C)(C)(C)") is None

    # test with invalid molecule
    assert dm.to_inchi_non_standard("") is None
    assert dm.to_inchikey_non_standard("") is None

    assert dm.to_inchi_non_standard("NC(N)=O") == "InChI=1/CH4N2O/c2-1(3)4/h(H4,2,3,4)/f/h2-3H2"
    assert dm.to_inchikey_non_standard("N=C(N)O") == "XSQUKJJJFZCRTK-ZIALIONUNA-N"


def test_to_smiles_with_indices():
    mol = dm.to_mol("Cn1c(=O)c2c(ncn2C)n(C)c1=O")

    # NOTE(hadim): not a great test I guess xD
    smiles = dm.to_smiles(mol, with_atom_indices=True)
    smiles = cast(str, smiles)
    assert all([str(i) in smiles for i in range(13)])
