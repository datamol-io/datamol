import pkg_resources

import pandas as pd


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
        data = pd.read_csv(f)
    return data
