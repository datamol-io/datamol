from typing import Union
from typing import Optional
from typing import List

import pathlib
import gzip

from rdkit import Chem
import pandas as pd
import fsspec

import datamol as dm


def read_csv(file_uri: Union[str, pathlib.Path], **kwargs):
    """Read a CSV file.

    NOTE(hadim): not sure this function add a real value. fsspec is also
    supported by pandas. Consider removing it.

    Args:
        file_uri (Union[str, pathlib.Path]): Path to the file. Can be a local
            path, HTTP, HTTPS, S3, GS, etc. See https://filesystem-spec.readthedocs.io/en/latest/
        kwargs: Arguments to pass to `pd.read_csv()`.

    Returns:
        df: a `pandas.DataFrame`
    """
    with fsspec.open(file_uri) as f:
        df = pd.read_csv(f, **kwargs)
    return df


def read_excel(
    file_uri: Union[str, pathlib.Path],
    sheet_name: Optional[Union[str, int, list]] = 0,
    **kwargs,
):
    """Read an excel file.

    NOTE(hadim): not sure this function add a real value. fsspec is also
    supported by pandas. Consider removing it.

    Args:
        file_uri (Union[str, pathlib.Path]): Path to the file. Can be a local
            path, HTTP, HTTPS, S3, GS, etc. See https://filesystem-spec.readthedocs.io/en/latest/
        sheet_name (Optional[Union[str, int, list]]): see `pandas.read_excel()` doc.
        kwargs: Arguments to pass to `pd.read_excel()`.

    Returns:
        df: a `pandas.DataFrame`
    """
    with fsspec.open(file_uri) as f:
        df = pd.read_excel(f, sheet_name=0, **kwargs)
    return df


def read_sdf(
    file_uri: Union[str, pathlib.Path],
    as_df: bool = False,
) -> Union[list, pd.DataFrame]:
    """Read an SDF file.

    Args:
        file_uri (Union[str, pathlib.Path]): Path to the file. Can be a local
            path, HTTP, HTTPS, S3, GS, etc. See https://filesystem-spec.readthedocs.io/en/latest/
        as_df (bool, optional): Whether to return a list mol or a pandas DataFrame. Default to False.

    Returns:
        df: a `pandas.DataFrame`
    """

    with fsspec.open(file_uri) as f:
        if str(file_uri).endswith(".gz") or str(file_uri).endswith(".gzip"):
            f = gzip.open(f)
        supplier = Chem.ForwardSDMolSupplier(f)
        mols = [mol for mol in supplier if mol is not None]

    if as_df:
        return dm.to_df(mols)

    return mols


def to_sdf(
    mols: Union[List[Chem.Mol], pd.DataFrame],
    file_uri: Union[str, pathlib.Path],
    smiles_column: Optional[Union[int, str]] = None,
):
    """Write molecules to a file.

    Args:
        mols (Union[List[Chem.Mol], pd.DataFrame]):
        file_uri (Union[str, pathlib.Path]): Path to the file. Can be a local
            path, HTTP, HTTPS, S3, GS, etc. See https://filesystem-spec.readthedocs.io/en/latest/
        smiles_column (Optional[Union[int, str]]): if `mols` is a dataframe, you must specify `smiles_column`
            for saving to sdf.

    Returns:
        df: a `pandas.DataFrame`
    """

    if isinstance(mols, pd.DataFrame):
        mols = dm.from_df(mols, smiles_column=smiles_column)

    with fsspec.open(file_uri, mode="w") as f:
        writer = Chem.SDWriter(f)
        for mol in mols:
            writer.write(mol)
        writer.close()
