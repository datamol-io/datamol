"""
The data module aims to provide a fast and convenient access to various molecular datasets.

---
"""

from typing import Optional
from typing import cast
from typing import Union
from typing import List
from typing import overload
from typing import Literal

import sys
import io
import functools

try:
    import importlib.resources as importlib_resources
except ImportError:
    import importlib_resources

import pandas as pd

from ..types import Mol
from ..io import read_sdf
from ..convert import from_df
from ..convert import render_mol_df


@functools.lru_cache()
def datamol_data_file_path(filename: str, dm_module: str = "datamol.data") -> str:
    if sys.version_info < (3, 9, 0):
        with importlib_resources.path(dm_module, filename) as p:
            data_path = p
    else:
        data_path = importlib_resources.files(dm_module).joinpath(filename)

    return str(data_path)


def open_datamol_data_file(
    filename: str,
    open_binary: bool = False,
    dm_module: str = "datamol.data",
):
    if sys.version_info < (3, 9, 0):
        if open_binary:
            file_context_manager = importlib_resources.open_binary(dm_module, filename)
        else:
            file_context_manager = importlib_resources.open_text(dm_module, filename)
    else:
        if open_binary:
            mode = "rb"
        else:
            mode = "r"

        file_context_manager = (
            importlib_resources.files(dm_module).joinpath(filename).open(mode=mode)
        )

    # NOTE(hadim): we assume the file always exists
    file_context_manager = cast(io.TextIOWrapper, file_context_manager)

    return file_context_manager


@overload
def freesolv(as_df: Literal[True] = True) -> pd.DataFrame: ...


@overload
def freesolv(as_df: Literal[False] = False) -> List[Mol]: ...


@overload
def freesolv(as_df: bool = True) -> Union[List[Mol], pd.DataFrame]: ...


def freesolv(as_df: bool = True) -> Union[List[Mol], pd.DataFrame]:
    """Return the FreeSolv dataset as a dataframe.

    The dataset contains 642 molecules and the following columns:
    `['iupac', 'smiles', 'expt', 'calc']`.

    Warning:
        This dataset is only meant to be used as a toy dataset for pedagogic and
        testing purposes. **It is not** a dataset for benchmarking, analysis or
        model training.
    """

    with open_datamol_data_file("freesolv.csv") as f:
        data = pd.read_csv(f)

    if not as_df:
        data = from_df(data)

    return data


@overload
def cdk2(as_df: Literal[True] = True, mol_column: Optional[str] = "mol") -> pd.DataFrame: ...


@overload
def cdk2(as_df: Literal[False] = False, mol_column: Optional[str] = "mol") -> List[Mol]: ...


@overload
def cdk2(
    as_df: bool = True, mol_column: Optional[str] = "mol"
) -> Union[List[Mol], pd.DataFrame]: ...


def cdk2(as_df: bool = True, mol_column: Optional[str] = "mol"):
    """Return the RDKit CDK2 dataset from `RDConfig.RDDocsDir, 'Book/data/cdk2.sdf'`.

    Args:
        as_df: Whether to return a list mol or a pandas DataFrame.
        mol_column: Name of the mol column. Only relevant if `as_df` is True.
    """

    with open_datamol_data_file("cdk2.sdf", open_binary=True) as f:
        data = read_sdf(f, as_df=as_df, mol_column=mol_column)
    return data


@overload
def solubility(as_df: Literal[True] = True, mol_column: Optional[str] = "mol") -> pd.DataFrame: ...


@overload
def solubility(as_df: Literal[False] = False, mol_column: Optional[str] = "mol") -> List[Mol]: ...


@overload
def solubility(
    as_df: bool = True, mol_column: Optional[str] = "mol"
) -> Union[List[Mol], pd.DataFrame]: ...


def solubility(as_df: bool = True, mol_column: Optional[str] = "mol"):
    """Return the RDKit solubility dataset from `RDConfig.RDDocsDir, 'Book/data/solubility.{train|test}.sdf'`.

    The dataframe or the list of molecules with contain a `split` column, either `train` or `test`.

    Args:
        as_df: Whether to return a list mol or a pandas DataFrame.
        mol_column: Name of the mol column. Only relevant if `as_df` is True.
    """

    with open_datamol_data_file("solubility.train.sdf", open_binary=True) as f:
        train = read_sdf(f, as_df=True, mol_column="mol", smiles_column=None)

    with open_datamol_data_file("solubility.test.sdf", open_binary=True) as f:
        test = read_sdf(f, as_df=True, mol_column="mol", smiles_column=None)

    train = cast(pd.DataFrame, train)
    test = cast(pd.DataFrame, test)

    train["split"] = "train"
    test["split"] = "test"

    # NOTE(hadim): LMAO RDkit consistency xD
    test = test.rename(columns={"SMILES": "smiles"})

    data = pd.concat([train, test], ignore_index=True)

    if as_df:
        if mol_column is None:
            data = data.drop(columns=["mol"])

        render_mol_df(data)
        return data

    return from_df(data, mol_column=mol_column)


@overload
def chembl_drugs(as_df: Literal[True] = True) -> pd.DataFrame: ...


@overload
def chembl_drugs(as_df: Literal[False] = False) -> List[Mol]: ...


def chembl_drugs(as_df: bool = True) -> Union[List[Mol], pd.DataFrame]:
    """A list of ~2.5k molecules from ChEMBL (all approved drugs) in SMILES format.
    Includes metadata indicating year of first approval, molecule chembl id, molecule type and pref_name.

    List was generated with ['Get_ChEMBL_Approved_Drugs.ipynb'](https://github.com/datamol-io/datamol/notebooks/Get_ChEMBL_Approved_Drugs.ipynb) on 2023-10-18.
    The notebook works with the chembl_webresource_client api to collect chembl IDs and metadata, then focuses on small molecules with valid SMILES and first approval date.
    """
    with open_datamol_data_file("chembl_approved_drugs.parquet", open_binary=True) as f:
        data = pd.read_parquet(f)

    if not as_df:
        data = from_df(data)

    return data


@overload
def chembl_samples(as_df: Literal[True] = True) -> pd.DataFrame: ...


@overload
def chembl_samples(as_df: Literal[False] = False) -> List[Mol]: ...


def chembl_samples(as_df: bool = True) -> Union[List[Mol], pd.DataFrame]:
    """A list of ~2k molecules from ChEMBL.

    Originally, proposed by Patrick Walters at <https://github.com/PatWalters/practical_cheminformatics_posts/tree/b4dae239a8b942dab3a41e637ac55d4491aee96f/molskill>.
    """

    with open_datamol_data_file("chembl_samples.csv") as f:
        data = pd.read_csv(f)

    if not as_df:
        data = from_df(data)

    return data
