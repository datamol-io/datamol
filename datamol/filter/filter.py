from typing import List

import pandas as pd

import datamol as dm

from .matches import n_matches
from .matches import get_descriptions
from .matches import get_props


def df_filter_smiles(
    list_of_smi: List[str],
    catalog_specifiers: List[str] = ["ALL"],
) -> pd.DataFrame:
    """See which smiles string has RDKit FilterCatalog hits

    Args:
        list_of_smi: A list of SMILES strings
        catalog_specifiers: Specify what kind of filtering catalog you want.

    Returns:
        df_filter_smiles: A `pd.DataFrame` table that contains the results of df_filter_smiles
    """
    des = []
    catalog = []
    num_hits = []
    scope = []
    for smi in list_of_smi:
        mol = dm.to_mol(smi)
        catalog.append(",".join(get_props(mol, catalog_specifiers=catalog_specifiers)))
        des.append(",".join(get_descriptions(mol, catalog_specifiers=catalog_specifiers)))
        num_hits.append(n_matches(mol, catalog_specifiers=catalog_specifiers))
        scope.append(",".join(get_props(mol, prop="Scope")))

    return pd.DataFrame(
        data=list(zip(list_of_smi, num_hits, catalog, des, scope)),
        columns=["SMILES", "num_of_matches", "Catalog", "Description", "Scope"],
    )


def filter_smiles(
    list_of_smi: List[str],
    catalog_specifiers: List[str] = ["ALL"],
) -> dict:
    """See which smiles string has RDKit FilterCatalog hits in a python dictionary

    Args:
        list_of_smi: A list of SMILES strings
        catalog_specifiers: Specify what kind of filtering catalog you want.

    Returns:
        filter_smiles: A python dictionary that contains the results of filter_smiles
    """
    des = []
    catalog = []
    num_hits = []
    for smi in list_of_smi:
        mol = dm.to_mol(smi)
        catalog.append(",".join(get_props(mol, catalog_specifiers=catalog_specifiers)))
        des.append(",".join(get_descriptions(mol, catalog_specifiers=catalog_specifiers)))
        num_hits.append(n_matches(mol, catalog_specifiers=catalog_specifiers))

    return dict(
        {"SMILES": list_of_smi, "num_of_matches": num_hits, "Catalog": catalog, "Description": des}
    )
