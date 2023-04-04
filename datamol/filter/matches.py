from typing import List

from rdkit.Chem import FilterCatalog


from .. import Mol

from .._version import is_lower_than_current_rdkit_version



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
    PAINS_params = FilterCatalog.FilterCatalogParams()
    PAINS_params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.names[filt_specifier])
    PAINS_catalog = FilterCatalog.FilterCatalog(PAINS_params) 

    entries = PAINS_catalog.GetMatches(mol)
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
    PAINS_params = FilterCatalog.FilterCatalogParams()
    PAINS_params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.names[filt_specifier])
    PAINS_catalog = FilterCatalog.FilterCatalog(PAINS_params) 

    entries = PAINS_catalog.GetMatches(mol)
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
    PAINS_params = FilterCatalog.FilterCatalogParams()
    PAINS_params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.names[filt_specifier])
    PAINS_catalog = FilterCatalog.FilterCatalog(PAINS_params) 

    entries = PAINS_catalog.GetMatches(mol)
    props = [entry.GetProp('FilterSet') for entry in entries]

    return props

