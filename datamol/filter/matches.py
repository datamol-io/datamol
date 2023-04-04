from typing import List

from rdkit.Chem import FilterCatalog


from .. import Mol

from .._version import is_lower_than_current_rdkit_version



def n_matches(    
    mol: Mol,
    filt_specifier: str = "ALL",
) -> int:
    """Compute the number of PAINS in a molecule.

    Rigid bonds are bonds that are not single and not in rings.

    Args:
        mol: A molecule.
        filter_specifier: Specify what kind of filtering catalog you want.

    Returns:
        n_matches_filter_catalog: number of matches based on specified catalog
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
    """Compute the number of PAINS in a molecule.

    Rigid bonds are bonds that are not single and not in rings.

    Args:
        mol: A molecule.
        filter_specifier: Specify what kind of filtering catalog you want.

    Returns:
        n_matches_filter_catalog: number of matches based on specified catalog
    """
    PAINS_params = FilterCatalog.FilterCatalogParams()
    PAINS_params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.names[filt_specifier])
    PAINS_catalog = FilterCatalog.FilterCatalog(PAINS_params) 

    entries = PAINS_catalog.GetMatches(mol)
    descriptions = [entry.GetDescription() for entry in entries]
    return descriptions

