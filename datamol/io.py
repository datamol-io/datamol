from typing import Union
from typing import Optional
from typing import List
from typing import Sequence

import io
import tempfile
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
        urlpath: Path to the file. Path to the file. Can be a local path, HTTP,
            HTTPS, S3, GS, etc. See https://filesystem-spec.readthedocs.io/en/latest/
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
        urlpath: Path to the file. Path to the file. Can be a local path, HTTP,
            HTTPS, S3, GS, etc. See https://filesystem-spec.readthedocs.io/en/latest/
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
        urlpath: Path to the file. Path to the file. Can be a local path, HTTP,
            HTTPS, S3, GS, etc. See https://filesystem-spec.readthedocs.io/en/latest/
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
        urlpath: Path to the file. Path to the file. Can be a local path, HTTP,
            HTTPS, S3, GS, etc. See https://filesystem-spec.readthedocs.io/en/latest/
        smiles_column: if `mols` is a dataframe, you must specify `smiles_column`
            for saving to sdf.
    """

    if isinstance(mols, pd.DataFrame):
        mols = dm.from_df(mols, smiles_column=smiles_column)

    # Filter out None values
    mols = [mol for mol in mols if mol is not None]

    with fsspec.open(urlpath, mode="w") as f:
        writer = Chem.SDWriter(f)
        for mol in mols:
            writer.write(mol)
        writer.close()


def to_smi(
    mols: Sequence[Chem.Mol],
    urlpath: Union[str, pathlib.Path],
    error_if_empty: bool = False,
):
    """Save a list of molecules in an `.smi` file.

    Args:
        mols: a list of molecules.
        urlpath: Path to the file. Path to the file. Can be a local path, HTTP,
            HTTPS, S3, GS, etc. See https://filesystem-spec.readthedocs.io/en/latest/
        error_if_empty: whether to raise and error if the input list is empty.
    """

    if len(mols) == 0 and error_if_empty:
        raise ValueError("The list of mols/smiles provided is empty.")

    # Filter out None values
    mols = [mol for mol in mols if mol is not None]

    with fsspec.open(urlpath, "w") as f:
        writer = Chem.SmilesWriter(f, includeHeader=False, nameHeader="")
        for mol in mols:
            writer.write(mol)
        writer.close()


def read_smi(
    urlpath: Union[str, pathlib.Path],
) -> Sequence[Chem.Mol]:
    """Read a list of smiles from am `.smi` file.

    Args:
        urlpath: Path to the file. Path to the file. Can be a local path, HTTP,
            HTTPS, S3, GS, etc. See https://filesystem-spec.readthedocs.io/en/latest/
    """

    active_path = urlpath

    # NOTE(hadim): the temporary local file copy
    # is because `SmilesMolSupplier` does not support
    # using file-like object, only path.

    # Copy to a local temporary path if the path is a remote one.
    if not fsspec.utils.can_be_local(str(urlpath)):
        active_path = pathlib.Path(tempfile.mkstemp()[1])
        dm.utils.copy_files(urlpath, active_path)

    # Read the molecules
    supplier = Chem.SmilesMolSupplier(str(active_path), titleLine=0)
    mols = [mol for mol in supplier if mol is not None]

    # Delete the local temporary path
    if not fsspec.utils.can_be_local(str(urlpath)):
        active_path.unlink()

    return mols
