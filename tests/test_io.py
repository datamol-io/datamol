import pathlib
import tempfile
import pytest
import io

from rdkit import Chem

import datamol as dm

ROOT_DIR = pathlib.Path(__file__).parent.resolve()
DATA_DIR = ROOT_DIR / "data"


def get_freesolv_csv_path():
    return DATA_DIR / "freesolv.csv"


def get_freesolv_excel_path():
    return DATA_DIR / "freesolv.xlsx"


def get_tubb3_sdf_path():
    return DATA_DIR / "TUBB3-observations.sdf"


def get_tubb3_sdf_gzip_path():
    return DATA_DIR / "TUBB3-observations.sdf.gz"


def test_read_csv():
    data_path = get_freesolv_csv_path()
    df = dm.read_csv(data_path)
    assert df.shape == (642, 4)
    assert list(df.columns) == ["iupac", "smiles", "expt", "calc"]


def test_read_excel():
    data_path = get_freesolv_excel_path()
    df = dm.read_excel(data_path, engine="openpyxl")
    assert df.shape == (642, 4)
    assert list(df.columns) == ["iupac", "smiles", "expt", "calc"]


def test_read_sdf():

    # sdf
    data_path = get_tubb3_sdf_path()
    mols = dm.read_sdf(data_path)
    assert len(mols) == 10
    for mol in mols:
        assert isinstance(mol, Chem.Mol)

    # sdf gzipped
    data_path = get_tubb3_sdf_gzip_path()
    mols = dm.read_sdf(data_path)
    assert len(mols) == 10
    for mol in mols:
        assert isinstance(mol, Chem.Mol)


def test_read_sdf_as_df():

    # sdf
    data_path = get_tubb3_sdf_path()
    df = dm.read_sdf(data_path, as_df=True)
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

    # sdf gzipped
    data_path = get_tubb3_sdf_gzip_path()
    df = dm.read_sdf(data_path, as_df=True)
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


def test_to_sdf():
    data_path = get_tubb3_sdf_gzip_path()

    # dataframe based
    df = dm.read_sdf(data_path, as_df=True)

    temp_path = tempfile.mktemp()
    dm.to_sdf(df, temp_path, smiles_column="smiles")

    new_df = dm.read_sdf(temp_path, as_df=True)
    assert df.equals(new_df)

    # list of mols based
    mols = dm.read_sdf(data_path, as_df=False)

    temp_path = tempfile.mktemp()
    dm.to_sdf(mols, temp_path)

    new_mols = dm.read_sdf(temp_path, as_df=False)
    assert [dm.to_smiles(mol) for mol in mols] == [dm.to_smiles(mol) for mol in new_mols]


def test_to_from_text():

    temp_file = pathlib.Path(tempfile.mkstemp()[1])

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
