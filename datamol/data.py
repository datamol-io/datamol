from typing import Optional

import os

import pkg_resources

import pandas as pd

from rdkit.Chem import RDConfig

from .io import read_sdf


def freesolv():
    """Return the FreeSolv dataset as a dataframe.

    The dataset contains 642 molecules and the following columns:
    `['iupac', 'smiles', 'expt', 'calc']`.

    Warning:
        This dataset is only meant to be used as a toy dataset for pedagogic and
        testing purposes. **It is not** a dataset for benchmarking, analysis or
        model training.
    """

    with pkg_resources.resource_stream("datamol", "data/freesolv.csv") as f:
        data = pd.read_csv(f)  # type: ignore
    return data


def cdk2(as_df: bool = True, mol_column: Optional[str] = "mol"):
    """Return the RDKit CDK2 dataset from `RDConfig.RDDocsDir, 'Book/data/cdk2.sdf'`.

    Args:
        as_df: Whether to return a list mol or a pandas DataFrame.
        mol_column: Name of the mol column. Only relevant if `as_df` is True.
    """

    return read_sdf(
        os.path.join(RDConfig.RDDocsDir, "Book/data/cdk2.sdf"),
        as_df=as_df,
        mol_column=mol_column,
    )
