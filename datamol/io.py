from typing import Union
from typing import Optional
from typing import List
from typing import Sequence

import pathlib
import gzip

from rdkit import Chem
import pandas as pd
import fsspec

import datamol as dm


def read_csv(urlpath: Union[str, pathlib.Path], **kwargs):
    """Read a CSV file.

    NOTE(hadim): not sure this function add a real value. fsspec is also
    supported by pandas. Consider removing it.

    Args:
        urlpath: Path to the file. Can be a local
            path, HTTP, HTTPS, S3, GS, etc. See https://filesystem-spec.readthedocs.io/en/latest/
        kwargs: Arguments to pass to `pd.read_csv()`.

    Returns:
        df: a `pandas.DataFrame`
    """
    with fsspec.open(urlpath) as f:
        df = pd.read_csv(f, **kwargs)
    return df


def read_excel(
    urlpath: Union[str, pathlib.Path],
    sheet_name: Optional[Union[str, int, list]] = 0,
    **kwargs,
):
    """Read an excel file.

    NOTE(hadim): not sure this function add a real value. fsspec is also
    supported by pandas. Consider removing it.

    Args:
        urlpath: Path to the file. Can be a local
            path, HTTP, HTTPS, S3, GS, etc. See https://filesystem-spec.readthedocs.io/en/latest/
        sheet_name: see `pandas.read_excel()` doc.
        kwargs: Arguments to pass to `pd.read_excel()`.

    Returns:
        df: a `pandas.DataFrame`
    """
    with fsspec.open(urlpath) as f:
        df = pd.read_excel(f, sheet_name=sheet_name, **kwargs)
    return df


def read_sdf(
    urlpath: Union[str, pathlib.Path],
    as_df: bool = False,
) -> Union[list, pd.DataFrame]:
    """Read an SDF file.

    Args:
        urlpath (Union[str, pathlib.Path]): Path to the file. Can be a local
            path, HTTP, HTTPS, S3, GS, etc. See https://filesystem-spec.readthedocs.io/en/latest/
        as_df: Whether to return a list mol or a pandas DataFrame. Default to False.

    Returns:
        df: a `pandas.DataFrame`
    """

    with fsspec.open(urlpath) as f:
        if str(urlpath).endswith(".gz") or str(urlpath).endswith(".gzip"):
            f = gzip.open(f)
        supplier = Chem.ForwardSDMolSupplier(f)
        mols = [mol for mol in supplier if mol is not None]

    if as_df:
        return dm.to_df(mols)

    return mols


def to_sdf(
    mols: Union[Sequence[Chem.Mol], pd.DataFrame],
    urlpath: Union[str, pathlib.Path],
    smiles_column: Optional[Union[int, str]] = None,
):
    """Write molecules to a file.

    Args:
        mols:
        urlpath: Path to the file. Can be a local
            path, HTTP, HTTPS, S3, GS, etc. See https://filesystem-spec.readthedocs.io/en/latest/
        smiles_column: if `mols` is a dataframe, you must specify `smiles_column`
            for saving to sdf.
    """

    if isinstance(mols, pd.DataFrame):
        mols = dm.from_df(mols, smiles_column=smiles_column)

    with fsspec.open(urlpath, mode="w") as f:
        writer = Chem.SDWriter(f)
        for mol in mols:
            if mol is not None:
                writer.write(mol)
        writer.close()


def to_text(
    mols: Sequence[Union[Chem.Mol, str]],
    urlpath: Union[str, pathlib.Path],
    error_if_empty: bool = False,
):
    """Save a list of molecules or smiles in a text file. One smiles per line.

    Args:
        mols: can be either a molecule or a smiles. If the mol to smiles conversion must be customized
            you should call `[dm.to_smiles(m, **custom_kargs) for m in mols]` before `dm.to_text`.
        urlpath: Path to the file. Can be a local
            path, HTTP, HTTPS, S3, GS, etc. See https://filesystem-spec.readthedocs.io/en/latest/
        error_if_empty: whether to raise and error if the input list is empty.
    """

    if len(mols) == 0 and error_if_empty:
        raise ValueError("The list of mols/smiles provided is empty.")

    if len(mols) >0 and isinstance(mols[0], Chem.Mol):
        mols = [dm.to_smiles(m) for m in mols]

    with fsspec.open(urlpath, mode="w") as f:
        for smiles in mols:
            f.write(f"{smiles}\n")


def read_text(
    urlpath: Union[str, pathlib.Path],
    as_mol: bool = True,
):
    """Read a list of smiles from a text file. One smiles per line.

    Args:
        urlpath: Path to the file. Can be a local
            path, HTTP, HTTPS, S3, GS, etc. See https://filesystem-spec.readthedocs.io/en/latest/
        as_mol: whether to convert the smiles to a list of molecules. If you need to use custom arguments,
            to convert to mol, use `[dm.to_mol(smiles, **custom_args) for smiles in smiles_list`. When set
            to false, smiles are loaded as-is without any validity check.
    """

    with fsspec.open(urlpath, mode="r") as f:
        mols = []
        for line in f.readlines():
            mols.append(line.strip())

    if as_mol:
        mols = [dm.to_mol(mol) for mol in mols]

    return mols
