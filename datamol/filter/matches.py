from typing import List

from rdkit.Chem import FilterCatalog


from .. import Mol

from .._version import is_lower_than_current_rdkit_version


def set_filter_params(
    filt_specifier
) -> object:
    """Utility function that sets filter set parameters


    Args:
        filter_specifier: Specify what kind of filtering catalog you want.

    Returns:
        Catalog: set filter set
    """
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

    matches = Catalog.GetMatches(mol)
    return len(matches)

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

def get_props(    
    mol: Mol,
    filt_specifier: str = "ALL",
    prop: str = 'FilterSet'
) -> List[str]:
    """Get the values of the props RDKit FilterCatalog hits for a molecule

    Args:
        mol: A molecule.
        filter_specifier: Specify what kind of filtering catalog you want.

    Returns:
        get_props: list of strings that expresses the values of the RDKit FilterCatalog hits
    """
    Catalog = set_filter_params(filt_specifier)
    entries = Catalog.GetMatches(mol)
    vals = [entry.GetProp(prop) for entry in entries]

    return vals

def get_prop_list(
    mol: Mol,
    filt_specifier: str = "ALL",
    deconstruct: bool = True,
) -> List[str]:
    """Get the list of props of RDKit FilterCatalog hits for a molecule

    Args:
        mol: A molecule.
        filter_specifier: Specify what kind of filtering catalog you want.
        deconstruct: Whether if you want to deconstruct prop list into more readable strings

    Returns:
        get_props_list: list of all prop keys for each RDKit FilterCatalog hit for a molecule
    """
    Catalog = set_filter_params(filt_specifier)
    entries = Catalog.GetMatches(mol)
    prop_list = [entry.GetPropList() for entry in entries]
    if deconstruct:
        decons_list = []
        for prop in prop_list:
            one_hit=[]
            for single in prop:
                one_hit.append(single)
            decons_list.append(one_hit)
        return decons_list
    return prop_list




# #
# def run_filter_catalog(
#         list_smi: List[str],
#         filt_specifier: str = "ALL",
#         num_of_threads: int = 1,
        
# ) -> List[str]:
#     """Run the filter catalog here on a list of smiles to find hits.

#     Args:
#         mol: A molecule.
#         filter_specifier: Specify what kind of filtering catalog you want.
#         deconstruct: Whether if you want to deconstruct prop list into more readable strings

#     Returns:
#         get_props_list: list of all prop keys for each RDKit FilterCatalog hit for a molecule
#     """
#     Catalog = set_filter_params(filt_specifier)
#     list_of_results = FilterCatalog.RunFilterCatalog(
#         filterCatalog=Catalog,
#         smiles=list_smi,
#         numThreads=num_of_threads,
#         )
#     for results in list_of_results:
#         for line in results:
#             print(line.GetDescription())
#     return results
        
