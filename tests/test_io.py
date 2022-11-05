import pytest
import io

from rdkit import Chem

import datamol as dm
import numpy as np
import pandas as pd


def test_read_csv(datadir):
    data_path = datadir / "freesolv.csv"
    df = dm.read_csv(data_path)
    assert df.shape == (642, 4)
    assert list(df.columns) == ["iupac", "smiles", "expt", "calc"]


def test_read_csv_with_mol_col(datadir):
    data_path = datadir / "freesolv.csv"
    df = dm.read_csv(data_path, smiles_column="smiles", mol_column="mol")
    assert df.shape == (642, 5)
    assert set(df.columns) == {"iupac", "smiles", "expt", "calc", "mol"}
    assert isinstance(df.iloc[0]["mol"], Chem.rdchem.Mol)


def test_read_excel(datadir):
    data_path = datadir / "freesolv.xlsx"
    df = dm.read_excel(data_path, engine="openpyxl")
    assert df.shape == (642, 4)
    assert set(df.columns) == {"iupac", "smiles", "expt", "calc"}


def test_read_excel_with_mol_col(datadir):
    data_path = datadir / "freesolv.xlsx"
    df = dm.read_excel(data_path, smiles_column="smiles", mol_column="mol")
    assert df.shape == (642, 5)
    assert set(df.columns) == {"iupac", "smiles", "expt", "calc", "mol"}
    assert isinstance(df.iloc[0]["mol"], Chem.rdchem.Mol)


def test_read_sdf(datadir):

    data_path = datadir / "TUBB3-observations.sdf"
    mols = dm.read_sdf(data_path)
    assert len(mols) == 10
    for mol in mols:
        assert isinstance(mol, Chem.rdchem.Mol)


def test_read_sdf_gz(datadir):

    data_path = datadir / "TUBB3-observations.sdf.gz"
    mols = dm.read_sdf(data_path)
    assert len(mols) == 10
    for mol in mols:
        assert isinstance(mol, Chem.rdchem.Mol)


def test_read_sdf_as_df(datadir):

    data_path = datadir / "TUBB3-observations.sdf"
    df = dm.read_sdf(data_path, as_df=True)
    assert df.shape == (10, 12)
    assert set(df.columns) == {
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
    }


def test_read_sdf_as_df_parallel(datadir):
    data_path = datadir / "TUBB3-observations.sdf"
    df = dm.read_sdf(data_path, as_df=True, n_jobs=1)
    df2 = dm.read_sdf(data_path, as_df=True, n_jobs=-1)
    pd.testing.assert_frame_equal(df, df2)


def test_read_sdf_as_df_mol_col(datadir):

    data_path = datadir / "TUBB3-observations.sdf"
    df = dm.read_sdf(data_path, as_df=True, mol_column="mol")
    assert "mol" in df.columns
    assert isinstance(df.iloc[0]["mol"], Chem.rdchem.Mol)


def test_read_sdf_gz_as_df(datadir):

    data_path = datadir / "TUBB3-observations.sdf.gz"
    df = dm.read_sdf(data_path, as_df=True)

    assert df.shape == (10, 12)
    assert set(df.columns) == {
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
    }


def test_to_sdf(datadir, tmp_path):
    data_path = datadir / "TUBB3-observations.sdf.gz"

    df = dm.read_sdf(data_path, as_df=True)

    sdf_path = tmp_path / "mols.sdf"
    dm.to_sdf(df, sdf_path, smiles_column="smiles")

    new_df = dm.read_sdf(sdf_path, as_df=True)
    assert df.equals(new_df)


def test_to_sdf_mols(datadir, tmp_path):
    data_path = datadir / "TUBB3-observations.sdf.gz"

    mols = dm.read_sdf(data_path, as_df=False)

    sdf_path = tmp_path / "mols.sdf"
    dm.to_sdf(mols, sdf_path)

    new_mols = dm.read_sdf(sdf_path, as_df=False)
    assert [dm.to_smiles(mol) for mol in mols] == [dm.to_smiles(mol) for mol in new_mols]


def test_to_from_text(tmp_path):

    temp_file = tmp_path / "mols.smi"

    smiles_list = [
        "Cn1c(=S)ccc2nc[nH]c21",
        "Clc1n[nH]c2c1=[NH+]C(c1ccc[nH+]c1)C[NH+]=2",
        "Fc1ccsc1",
        "N#Cc1cc2c(o1)[NH2+]CCN2Cn1cnc2c1CSCC2",
        "O=CN1CCC2NC=CC2C1",
        "Oc1[nH]nc2c1-n1ncnc1C2",
        "OC1=NNC2(OC=CCO2)C2(C3CCCc4nonc43)NN=NN12",
        "[NH-]Sc1cc2nc[nH+]cc2o1",
        "[NH3+]C12CNCCOC1(N1CCCCC1)C=C(F)NC2",
    ]
    mols = [dm.to_mol(m) for m in smiles_list]

    # Save from text and read from text
    dm.to_smi(mols, temp_file)
    loaded_mols = dm.read_smi(temp_file)
    loaded_smiles = [dm.to_smiles(m) for m in loaded_mols]
    assert loaded_smiles == smiles_list

    # Check error raised when list is empty
    with pytest.raises(ValueError):
        dm.to_smi([], temp_file, error_if_empty=True)

    temp_file.unlink()

    # Check file like object works too
    file_like = io.StringIO()
    dm.to_smi(mols, file_like)
    assert file_like.getvalue().strip().split("\n") == smiles_list


def test_to_sdf_single_mol(tmp_path):

    sdf_path = tmp_path / "test.sdf"

    smiles = "CC1(C2C(C3C(C(=O)C(=C(C3(C(=O)C2=C(C4=C1C=CC=C4O)O)O)O)C(=O)N)N(C)C)O)O"
    mol = dm.to_mol(smiles)
    dm.to_sdf(mol, sdf_path)

    mols = dm.read_sdf(sdf_path)
    assert dm.to_smiles(mol) == dm.to_smiles(mols[0])


def test_sdf_props_and_conformer_preserved(tmp_path):

    sdf_path = tmp_path / "test.sdf"

    # Generate an SDF file
    props = dict(test_int=588, test_str="hello")
    smiles = "CC1(C2C(C3C(C(=O)C(=C(C3(C(=O)C2=C(C4=C1C=CC=C4O)O)O)O)C(=O)N)N(C)C)O)O"

    mol = dm.to_mol(smiles)
    mol = dm.set_mol_props(mol, props)
    mol = dm.conformers.generate(mol, n_confs=1)
    pos = mol.GetConformer().GetPositions()
    dm.to_sdf(mol, sdf_path)

    # Read sdf file
    mols = dm.read_sdf(sdf_path)
    mol = mols[0]

    # Check properties
    assert mol.GetPropsAsDict() == props

    # Check conformer
    conf = mol.GetConformer()
    assert mol.GetNumConformers() == 1
    assert conf.Is3D()
    np.testing.assert_almost_equal(conf.GetPositions(), pos, decimal=4)


def test_read_save_molblock():
    mol = dm.to_mol("Cn1c(=O)c2c(ncn2C)n(C)c1=O")

    # to molblock
    molblock = dm.to_molblock(mol)

    assert isinstance(molblock, str)
    assert "END" in molblock
    assert "V2000" in molblock
    assert "RDKit" in molblock

    # read molblock
    mol2 = dm.read_molblock(molblock)
    assert dm.same_mol(mol, mol2)


def test_read_molblock_invalid():

    mol = dm.read_molblock("hello")
    assert mol is None

    with pytest.raises(ValueError):
        dm.read_molblock("hello", fail_if_invalid=True)


def test_to_xlsx(tmp_path):
    excel_path1 = tmp_path / "test1.xlsx"
    excel_path2 = tmp_path / "test2.xlsx"

    data = dm.freesolv()
    data = data.iloc[:10]
    data["mol"] = data["smiles"].apply(dm.to_mol)

    # write from df
    dm.to_xlsx(data, excel_path1)
    assert excel_path1.exists()

    # write from list of molecules
    mols = dm.from_df(data)
    dm.to_xlsx(mols, excel_path2)
    assert excel_path2.exists()


def test_to_xlsx_empty():
    mols = [None]
    with pytest.raises(ValueError):
        dm.to_xlsx(mols, "/dev/null")  # type: ignore
