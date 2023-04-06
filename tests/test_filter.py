import pytest

from rdkit.Chem import FilterCatalog

import pandas as pd
import datamol as dm


list_of_smi = ['O=C(Cn1cnc2c1c(=O)n(C)c(=O)n2C)N/N=C/c1c(O)ccc2c1cccc2','O=C(C)Oc1ccccc1C(=O)O','CCCCC']

def test_set_filter_params():

    assert isinstance(dm.filter.matches.set_filter_params(['ALL']), FilterCatalog.FilterCatalog)

    #cases where there are duplicates in your list of specfied catalogs
    with pytest.raises(ValueError):
        dm.filter.matches.set_filter_params(['PAINS', 
                                             'PAINS'])

    #a case where there are duplicates but the word 'PAINS' still exists
    with pytest.raises(ValueError):
        dm.filter.matches.set_filter_params(['PAINS_A',
                                             'PAINS_A',
                                             'PAINS'])

    #A case where everythign is unique but 'PAINS' is in the list 
    #with other pains
    Catalog = dm.filter.matches.set_filter_params(['PAINS_A','PAINS','PAINS_B'])
    assert Catalog.GetNumEntries == 817
