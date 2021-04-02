# type: ignore
import pathlib
import unittest

import pytest

import pandas as pd
from rdkit import Chem

import datamol as dm

ROOT_DIR = pathlib.Path(__file__).parent.resolve()
DATA_DIR = ROOT_DIR / "data"


class TestConvert(unittest.TestCase):
    def get_freesolv_csv_path(self):
        return DATA_DIR / "freesolv.csv"

    def get_freesolv_excel_path(self):
        return DATA_DIR / "freesolv.xlsx"

    def get_tubb3_sdf_path(self):
        return DATA_DIR / "TUBB3-observations.sdf"

    def get_tubb3_sdf_gzip_path(self):
        return DATA_DIR / "TUBB3-observations.sdf.gz"

    def test_to_smiles(self):

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

    def test_to_selfies(self):
        smiles = "CC(=O)Oc1ccccc1C(=O)O"
        mol = dm.to_mol(smiles)

        true_sf = "[C][C][Branch1_2][C][=O][O][C][=C][C][=C][C][=C][Ring1][Branch1_2][C][Branch1_2][C][=O][O]"

        selfies = dm.to_selfies(smiles)
        assert selfies == true_sf

        selfies = dm.to_selfies(mol)
        assert selfies == true_sf

    def test_from_selfies(self):
        selfies = "[C][C][Branch1_2][C][=O][O][C][=C][C][=C][C][=C][Ring1][Branch1_2][C][Branch1_2][C][=O][O]"

        smiles = dm.from_selfies(selfies, as_mol=False)
        assert smiles == "CC(=O)OC1=CC=CC=C1C(=O)O"

        mol = dm.from_selfies(selfies, as_mol=True)
        assert dm.to_smiles(mol) == "CC(=O)Oc1ccccc1C(=O)O"

    def test_to_smarts(self):
        smiles = "O=C(C)Oc1ccccc1C(=O)O"
        mol = dm.to_mol(smiles)

        smarts = dm.to_smarts(mol, keep_hs=True)
        assert smarts == "[CH3]-[C](=[O])-[O]-[c]1:[cH]:[cH]:[cH]:[cH]:[c]:1-[C](=[O])-[OH]"

        smarts = dm.to_smarts(mol, keep_hs=False)
        assert smarts == "[CH3]-[C](=[O])-[O]-[c]1:[cH]:[cH]:[cH]:[cH]:[c]:1-[C](=[O])-[OH]"

        assert dm.to_smarts(None) is None

    def test_inchi(self):
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

    def test_to_df(self):
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

    def test_to_df(self):

        data_path = self.get_tubb3_sdf_path()
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

    def test_from_df(self):
        data_path = self.get_tubb3_sdf_path()
        df = dm.read_sdf(data_path, as_df=True)

        mols = dm.from_df(df)

        assert len(mols) == 10
        assert isinstance(mols[0], Chem.Mol)

        assert set(mols[0].GetPropsAsDict().keys()) == set(
            [
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
        )

        assert dm.from_df(pd.DataFrame()) == []

    def test_to_cxsmiles(self):
        mol = dm.to_mol("OC1=CC2CCCCC2[N:1]=C1")
        smiles = dm.to_smiles(mol, cxsmiles=True)
        assert smiles == "OC1=CC2CCCCC2[N:1]=C1 |atomProp:9.molAtomMapNumber.1|"

    def test_to_smiles_fail(self):
        smiles = dm.to_smiles(55, allow_to_fail=False)
        assert smiles == None

        # NOTE(hadim): ideally you want to catch only `Boost.Python.ArgumentError` here.
        with pytest.raises(Exception):
            dm.to_smiles(55, allow_to_fail=True)
