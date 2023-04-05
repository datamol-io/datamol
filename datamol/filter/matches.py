from typing import List

from rdkit.Chem import FilterCatalog


from .. import Mol

from .._version import is_lower_than_current_rdkit_version


def set_filter_params(
    filt_specifier
) -> object:
    Catalog = FilterCatalog.FilterCatalogParams()
    Catalog.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.names[filt_specifier])
    Catalog = FilterCatalog.FilterCatalog(Catalog) 
    return Catalog


def n_matches(    
    mol: Mol,
    filt_specifier: str = "ALL",
) -> int:
    """Compute the number of RDKIT FilterCatalogs hits in a molecule.


    Args:
        mol: A molecule.
        filter_specifier: Specify what kind of filtering catalog you want.

    Returns:
        n_matches: number of matches based on specified catalog
    """
    Catalog = set_filter_params(filt_specifier)

    entries = Catalog.GetMatches(mol)
    return len(entries)

def get_descriptions(    
    mol: Mol,
    filt_specifier: str = "ALL",
) -> List[str]:
    """Get the descriptions of RDKIT FilterCatalogs hits in a molecule.


    Args:
        mol: A molecule.
        filter_specifier: Specify what kind of filtering catalog you want.

    Returns:
        get_descriptions: list of strings that expresses the description of each hit
    """
    Catalog = set_filter_params(filt_specifier)

    entries = Catalog.GetMatches(mol)
    descriptions = [entry.GetDescription() for entry in entries]
    return descriptions

def get_filter_set(    
    mol: Mol,
    filt_specifier: str = "ALL",
) -> List[str]:
    """Get the filter set of the RDKIT FilterCatalogs hits in a molecule


    Args:
        mol: A molecule.
        filter_specifier: Specify what kind of filtering catalog you want.

    Returns:
        get_descriptions: list of strings that expresses the filter set of each hit
    """
    Catalog = set_filter_params(filt_specifier)

    entries = Catalog.GetMatches(mol)
    props = [entry.GetProp('FilterSet') for entry in entries]

    return props

