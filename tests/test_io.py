import pytest
import io
import gzip

from rdkit import Chem

import datamol as dm
import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal


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

    # Read all sdf
    mols = dm.read_sdf(data_path)
    assert len(mols) == 10
    for mol in mols:
        assert isinstance(mol, Chem.rdchem.Mol)

    # Read 5 molecules
    mols = dm.read_sdf(data_path, max_num_mols=5)
    assert len(mols) == 5
    for mol in mols:
        assert isinstance(mol, Chem.rdchem.Mol)

    # Read more than the max number of mols in the file
    mols = dm.read_sdf(data_path, max_num_mols=111)
    assert len(mols) == 10
    for mol in mols:
        assert isinstance(mol, Chem.rdchem.Mol)

    data_path = datadir / "TUBB3-observations-last-broken.sdf"

    # Read all sdf with last mol being broken
    mols = dm.read_sdf(data_path)
    assert len(mols) == 9
    for mol in mols:
        assert isinstance(mol, Chem.rdchem.Mol)

    # Read all sdf with last mol being broken
    mols = dm.read_sdf(data_path, discard_invalid=False)
    assert len(mols) == 10
    for mol in mols[:-1]:
        assert isinstance(mol, Chem.rdchem.Mol)
    assert mols[-1] is None

    # Read all sdf with last mol being broken
    df = dm.read_sdf(data_path, discard_invalid=False, as_df=True, mol_column="mols")
    assert len(mols) == 10
    for mol in df["mols"].iloc[:-1]:
        assert isinstance(mol, Chem.rdchem.Mol)
    assert df["mols"].iloc[-1] is None


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


def test_read_mol2(datadir):
    data_path = datadir / "test.mol2"

    # to list of mols
    mols = dm.read_mol2file(data_path)

    assert isinstance(mols[0], Chem.rdchem.Mol)
    assert isinstance(mols[1], Chem.rdchem.Mol)
    assert isinstance(mols[2], Chem.rdchem.Mol)
    # cases where mol2 formats are damaged
    assert mols[3] is None
    assert mols[4] is None
    assert mols[5] is None
    assert mols[6] is None
    assert mols[7] is None

    firstMol = dm.to_mol("c1ccncc1")
    secondMol = dm.to_mol("c1c[nH]cn1")

    assert dm.same_mol(mols[0], firstMol)
    assert dm.same_mol(mols[1], secondMol)
    assert dm.same_mol(mols[2], secondMol)

    # a case where exception is raised because of None values
    with pytest.raises(ValueError):
        mols = dm.read_mol2file(data_path, fail_if_invalid=True)


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


def test_read_pdbblock():
    pdbblock = """HETATM    1  N1  UNL     1       3.707  -1.649   1.733  1.00  0.00           N
HETATM    2  C1  UNL     1       2.987  -1.188   0.938  1.00  0.00           C
HETATM    3  C2  UNL     1       1.869  -0.581   0.250  1.00  0.00           C
HETATM    4  C3  UNL     1       1.632   0.783   0.355  1.00  0.00           C
HETATM    5  C4  UNL     1       0.344   1.270   0.178  1.00  0.00           C
HETATM    6  C5  UNL     1      -0.755   0.429   0.049  1.00  0.00           C
HETATM    7  N2  UNL     1      -1.944   0.727   0.745  1.00  0.00           N
HETATM    8  C6  UNL     1      -3.233   0.937   0.173  1.00  0.00           C
HETATM    9  C7  UNL     1      -4.177   1.582   1.134  1.00  0.00           C
HETATM   10  C8  UNL     1      -5.125   2.560   0.500  1.00  0.00           C
HETATM   11  C9  UNL     1      -4.297   3.560  -0.257  1.00  0.00           C
HETATM   12  C10 UNL     1      -3.440   2.952  -1.324  1.00  0.00           C
HETATM   13  C11 UNL     1      -3.253   1.476  -1.208  1.00  0.00           C
HETATM   14  N3  UNL     1      -0.470  -0.872  -0.207  1.00  0.00           N
HETATM   15  C12 UNL     1       0.804  -1.335  -0.226  1.00  0.00           C
HETATM   16  N4  UNL     1       0.815  -2.626  -0.634  1.00  0.00           N
HETATM   17  C13 UNL     1      -0.473  -3.035  -0.715  1.00  0.00           C
HETATM   18  C14 UNL     1      -0.971  -4.331  -0.663  1.00  0.00           C
HETATM   19  C15 UNL     1      -2.225  -4.591  -0.133  1.00  0.00           C
HETATM   20  C16 UNL     1      -3.085  -3.511   0.005  1.00  0.00           C
HETATM   21  C17 UNL     1      -2.511  -2.281   0.175  1.00  0.00           C
HETATM   22  C18 UNL     1      -1.271  -1.954  -0.359  1.00  0.00           C
HETATM   23  C19 UNL     1       2.789   1.645   0.539  1.00  0.00           C
HETATM   24  C20 UNL     1       3.510   2.001  -0.595  1.00  0.00           C
HETATM   25  C21 UNL     1       4.873   2.200  -0.578  1.00  0.00           C
HETATM   26  C22 UNL     1       5.560   2.162   0.605  1.00  0.00           C
HETATM   27  C23 UNL     1       4.838   1.985   1.767  1.00  0.00           C
HETATM   28  C24 UNL     1       3.502   1.683   1.730  1.00  0.00           C
CONECT    1    2    2    2
CONECT    2    3
CONECT    3    4    4   15
CONECT    4    5   23
CONECT    5    6    6
CONECT    6    7   14
CONECT    7    8
CONECT    8    9   13
CONECT    9   10
CONECT   10   11
CONECT   11   12
CONECT   12   13
CONECT   14   15   22
CONECT   15   16   16
CONECT   16   17
CONECT   17   18   18   22
CONECT   18   19
CONECT   19   20   20
CONECT   20   21
CONECT   21   22   22
CONECT   23   24   24   28
CONECT   24   25
CONECT   25   26   26
CONECT   26   27
CONECT   27   28   28
END"""

    mol = dm.read_pdbblock(pdbblock)

    print(dm.to_smiles(mol))

    assert mol is not None
    assert mol.GetNumAtoms() == 28
    assert dm.to_inchikey(mol) == "ZVAMKEUGOUZEJZ-UHFFFAOYSA-N"

    conf = mol.GetConformer()
    assert conf.Is3D()


def test_to_pdbblock():
    mol = dm.to_mol("N#Cc1c(cc(NC2CCCCC2)n2c1nc1ccccc21)-c1ccccc1")
    molblock = dm.to_pdbblock(mol)

    assert "HETATM" in molblock
    assert "CONECT" in molblock
    assert "C19" in molblock
    assert "C11" in molblock

    mol2 = dm.read_pdbblock(molblock)
    assert dm.to_inchikey(mol2) == dm.to_inchikey(mol)


def test_read_pdbfile(tmp_path):
    mol = dm.to_mol("N#Cc1c(cc(NC2CCCCC2)n2c1nc1ccccc21)-c1ccccc1")
    mol = dm.conformers.generate(mol, n_confs=1)

    pdb_path = tmp_path / "test.pdb"
    dm.to_pdbfile(mol, pdb_path)

    mol2 = dm.read_pdbfile(pdb_path)
    assert dm.to_inchikey(mol2) == dm.to_inchikey(mol)


def test_to_pdbfile(tmp_path):
    mol = dm.to_mol("N#Cc1c(cc(NC2CCCCC2)n2c1nc1ccccc21)-c1ccccc1")
    mol = dm.conformers.generate(mol, n_confs=1)

    pdb_path = tmp_path / "test.pdb"

    dm.to_pdbfile(mol, pdb_path)

    with open(pdb_path) as f:
        molblock = f.read()

    assert "HETATM" in molblock
    assert "CONECT" in molblock
    assert "C19" in molblock
    assert "C11" in molblock


# tests begining open_df and save_df
@pytest.fixture
def tmp_pathIO(tmpdir):
    return str(tmpdir.mkdir("temp").join("test_file"))


def test_dataframe_csv(tmp_pathIO):
    data = pd.DataFrame({"col1": [1, 2, 3], "col2": ["a", "b", "c"]})
    dm.save_df(data, tmp_pathIO + ".csv")
    loaded_data = dm.open_df(tmp_pathIO + ".csv")
    assert_frame_equal(data, loaded_data)


def test_dataframe_excel(tmp_pathIO):
    data = pd.DataFrame({"col1": [1, 2, 3], "col2": ["a", "b", "c"]})
    dm.save_df(data, tmp_pathIO + ".xlsx")
    loaded_data = dm.open_df(tmp_pathIO + ".xlsx")
    assert_frame_equal(data, loaded_data)


def test_dataframe_parquet(tmp_pathIO):
    data = pd.DataFrame({"col1": [1, 2, 3], "col2": ["a", "b", "c"]})
    dm.save_df(data, tmp_pathIO + ".parquet")
    loaded_data = dm.open_df(tmp_pathIO + ".parquet")
    assert_frame_equal(data, loaded_data)


def test_dataframe_json(tmp_pathIO):
    data = pd.DataFrame({"col1": [1, 2, 3], "col2": ["a", "b", "c"]})
    dm.save_df(data, tmp_pathIO + ".json")
    loaded_data = dm.open_df(tmp_pathIO + ".json")
    assert_frame_equal(data, loaded_data)


def test_dataframe_sdf(tmp_pathIO):
    data = pd.DataFrame({"smiles": ["CC", "CCCC"], "col2": ["a", "b"]})
    dm.save_df(data, tmp_pathIO + ".sdf")
    loaded_data = dm.open_df(tmp_pathIO + ".sdf")
    assert_frame_equal(data, loaded_data)


def test_dataframe_sdf_gz(tmp_pathIO):
    data = pd.DataFrame({"smiles": ["CC", "CCCC"], "col2": ["a", "b"]})
    dm.save_df(data, tmp_pathIO + ".sdf")
    with open(tmp_pathIO + ".sdf", "rb") as f_in:
        with gzip.open(tmp_pathIO + ".sdf.gz", "wb") as f_out:
            f_out.writelines(f_in)
    loaded_data = dm.open_df(tmp_pathIO + ".sdf.gz")
    assert_frame_equal(data, loaded_data)


def test_save_dataframe_invalid(tmp_pathIO):
    data = pd.DataFrame({"col1": [1, 2, 3], "col2": ["a", "b", "c"]})
    with pytest.raises(ValueError):
        dm.save_df(data, tmp_pathIO + ".invalid")


def test_load_dataframe_invalid(tmp_pathIO):
    with pytest.raises(ValueError):
        dm.open_df(tmp_pathIO + ".invalid")


# No assertions, just checking that the function runs with kwargs
def test_file_with_kwargs(tmp_pathIO):
    data = pd.DataFrame({"col1": [1, 2, 3], "col2": ["a", "b", "c"]})
    dm.save_df(data, tmp_pathIO + ".csv", doublequote=False, sep=";")
    dm.open_df(tmp_pathIO + ".csv", verbose=True)
