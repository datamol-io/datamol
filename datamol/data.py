import pkg_resources

import pandas as pd


def freesolv():
    with pkg_resources.resource_stream("datamol", "data/freesolv.csv") as f:
        data = pd.read_csv(f)
    return data
