import unittest
import pathlib
import tempfile

from rdkit import Chem

import datamol as dm

ROOT_DIR = pathlib.Path(__file__).parent.resolve()
DATA_DIR = ROOT_DIR / "data"


class TestIO(unittest.TestCase):
    def get_freesolv_csv_path(self):
        return DATA_DIR / "freesolv.csv"

    def get_freesolv_excel_path(self):
        return DATA_DIR / "freesolv.xlsx"

    def get_tubb3_sdf_path(self):
        return DATA_DIR / "TUBB3-observations.sdf"

    def get_tubb3_sdf_gzip_path(self):
        return DATA_DIR / "TUBB3-observations.sdf.gz"

    def test_read_csv(self):
        data_path = self.get_freesolv_csv_path()
        df = dm.read_csv(data_path)
        assert df.shape == (642, 4)
        assert list(df.columns) == ["iupac", "smiles", "expt", "calc"]

    def test_read_excel(self):
        data_path = self.get_freesolv_excel_path()
        df = dm.read_excel(data_path, engine="openpyxl")
        assert df.shape == (642, 4)
        assert list(df.columns) == ["iupac", "smiles", "expt", "calc"]

    def test_read_sdf(self):

        # sdf
        data_path = self.get_tubb3_sdf_path()
        mols = dm.read_sdf(data_path)
        assert len(mols) == 10
        for mol in mols:
            assert isinstance(mol, Chem.Mol)

        # sdf gzipped
        data_path = self.get_tubb3_sdf_gzip_path()
        mols = dm.read_sdf(data_path)
        assert len(mols) == 10
        for mol in mols:
            assert isinstance(mol, Chem.Mol)

    def test_read_sdf_as_df(self):

        # sdf
        data_path = self.get_tubb3_sdf_path()
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
        data_path = self.get_tubb3_sdf_gzip_path()
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

    def test_to_sdf(self):
        data_path = self.get_tubb3_sdf_gzip_path()

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
