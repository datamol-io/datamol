import pytest

from rdkit.Chem import FilterCatalog
import pandas as pd
import datamol as dm

list_of_smi = [
    "O=C(Cn1cnc2c1c(=O)n(C)c(=O)n2C)N/N=C/c1c(O)ccc2c1cccc2",
    "O=C(C)Oc1ccccc1C(=O)O",
    "CCCCC",
]
n_of_matches = [3, 1, 0]
l_of_props = [["PAINS_A", "Brenk", "Brenk"], ["Brenk"], []]
l_of_descriptions = [
    ["hzone_phenol_A(479)", "imine_1", "Oxygen-nitrogen_single_bond"],
    ["phenol_ester"],
    [],
]
l_of_prop_lists = [
    [
        ["description", "FilterSet", "Reference", "Scope"],
        ["description", "FilterSet", "Reference", "Scope"],
        ["description", "FilterSet", "Reference", "Scope"],
    ],
    [["description", "FilterSet", "Reference", "Scope"]],
    [],
]
# set number of entries in each RDKIT filter sets.
PAINS_A = 16
PAINS_B = 55
PAINS_C = 409
PAINS = PAINS_A + PAINS_B + PAINS_C
BRENK = 105
ZINC = 50
NIH = 180
ALL = PAINS + BRENK + ZINC + NIH


def test_set_filter_params():
    assert isinstance(dm.filter.matches.set_filter_params(["ALL"]), FilterCatalog.FilterCatalog)

    # Here we are trying to test the set_filter_params handeling of
    # the user-defined list that specifies filter sets.
    Catalog = dm.filter.matches.set_filter_params(["PAINS_A"])
    assert Catalog.GetNumEntries() == PAINS_A
    Catalog = dm.filter.matches.set_filter_params(["PAINS_B"])
    assert Catalog.GetNumEntries() == PAINS_B
    Catalog = dm.filter.matches.set_filter_params(["PAINS_C"])
    assert Catalog.GetNumEntries() == PAINS_C
    Catalog = dm.filter.matches.set_filter_params(["PAINS"])
    assert Catalog.GetNumEntries() == PAINS
    Catalog = dm.filter.matches.set_filter_params(["BRENK"])
    assert Catalog.GetNumEntries() == BRENK
    Catalog = dm.filter.matches.set_filter_params(["ZINC"])
    assert Catalog.GetNumEntries() == ZINC
    Catalog = dm.filter.matches.set_filter_params(["NIH"])
    assert Catalog.GetNumEntries() == NIH
    Catalog = dm.filter.matches.set_filter_params(["ALL"])
    assert Catalog.GetNumEntries() == ALL

    # cases where there are duplicates in your list of specfied catalogs
    Catalog = dm.filter.matches.set_filter_params(["PAINS", "PAINS"])
    assert Catalog.GetNumEntries() == PAINS

    # a case where there are duplicates but the word 'PAINS' still exists
    Catalog = dm.filter.matches.set_filter_params(["PAINS_A", "PAINS_A", "PAINS"])
    assert Catalog.GetNumEntries() == PAINS
    # A case where everythign is unique but 'PAINS' is in the list
    # with other pains
    Catalog = dm.filter.matches.set_filter_params(["PAINS_A", "PAINS", "PAINS_B"])
    assert Catalog.GetNumEntries() == PAINS

    # Where all three types of PAINS family are in one list
    Catalog = dm.filter.matches.set_filter_params(["PAINS_A", "PAINS_B", "PAINS_C"])
    assert Catalog.GetNumEntries() == PAINS

    # Where all three types of PAINS family are in one list with duplicates
    Catalog = dm.filter.matches.set_filter_params(["PAINS_A", "PAINS_A", "PAINS_B", "PAINS_C"])
    assert Catalog.GetNumEntries() == PAINS

    Catalog = dm.filter.matches.set_filter_params(
        ["PAINS_A", "PAINS_A", "PAINS_B", "PAINS_B", "PAINS_C", "PAINS_C"]
    )
    assert Catalog.GetNumEntries() == PAINS

    # Then different combinations of others.
    Catalog = dm.filter.matches.set_filter_params(["PAINS", "NIH"])
    assert Catalog.GetNumEntries() == PAINS + NIH

    Catalog = dm.filter.matches.set_filter_params(["PAINS", "NIH", "BRENK"])
    assert Catalog.GetNumEntries() == PAINS + NIH + BRENK

    Catalog = dm.filter.matches.set_filter_params(["PAINS", "NIH", "BRENK", "ZINC"])
    assert Catalog.GetNumEntries() == PAINS + NIH + BRENK + ZINC

    Catalog = dm.filter.matches.set_filter_params(["PAINS_A", "NIH", "BRENK", "ZINC"])
    assert Catalog.GetNumEntries() == PAINS_A + NIH + BRENK + ZINC

    Catalog = dm.filter.matches.set_filter_params(["PAINS", "NIH", "BRENK", "ZINC", "ALL"])
    assert Catalog.GetNumEntries() == ALL

    # negative cases
    with pytest.raises(ValueError):
        dm.filter.matches.set_filter_params([])
    with pytest.raises(ValueError):
        dm.filter.matches.set_filter_params([""])
    with pytest.raises(ValueError):
        dm.filter.matches.set_filter_params(
            [
                "",
            ]
        )
    with pytest.raises(ValueError):
        dm.filter.matches.set_filter_params(["", ""])

    # If the set strings have lower case
    Catalog = dm.filter.matches.set_filter_params(["pains_a", "paiNS_b"])
    assert Catalog.GetNumEntries() == PAINS_A + PAINS_B

    with pytest.raises(ValueError):
        dm.filter.matches.set_filter_params(["pains_a", ""])

    with pytest.raises(KeyError):
        dm.filter.matches.set_filter_params(["pain", "pains_b"])

    with pytest.raises(ValueError):
        dm.filter.matches.set_filter_params(["pain", ""])


def test_n_matches():
    for i, smi in enumerate(list_of_smi, 0):
        mol = dm.to_mol(smi)
        assert dm.filter.n_matches(mol) == n_of_matches[i]


def test_get_description():
    for i, smi in enumerate(list_of_smi, 0):
        mol = dm.to_mol(smi)
        assert dm.filter.get_descriptions(mol) == l_of_descriptions[i]


def test_get_props():
    for i, smi in enumerate(list_of_smi, 0):
        mol = dm.to_mol(smi)
        assert dm.filter.get_props(mol) == l_of_props[i]


def test_get_prop_list():
    for i, smi in enumerate(list_of_smi, 0):
        mol = dm.to_mol(smi)
        assert dm.filter.get_prop_list(mol) == l_of_prop_lists[i]


def test_run_filter_catalog():
    l_of_results = dm.filter.run_filter_catalog(list_of_smi)
    for i, result in enumerate(l_of_results, 0):
        for j, res in enumerate(result, 0):
            assert res.GetDescription() == l_of_descriptions[i][j]
    assert isinstance(l_of_results, FilterCatalog.FilterCatalogListOfEntryList)


def test_df_filter_smiles():
    df = pd.DataFrame()
    df = dm.filter.df_filter_smiles(list_of_smi)

    scope = []
    catalog = []
    des = []
    for i, smi in enumerate(list_of_smi, 0):
        mol = dm.to_mol(smi)
        scope.append(",".join(dm.filter.get_props(mol, prop="Scope")))
        catalog.append(",".join(l_of_props[i]))
        des.append(",".join(l_of_descriptions[i]))
    df_t = pd.DataFrame(
        data=list(zip(list_of_smi, n_of_matches, catalog, des, scope)),
        columns=["SMILES", "num_of_matches", "Catalog", "Description", "Scope"],
    )

    assert isinstance(df, pd.DataFrame)
    assert isinstance(df_t, pd.DataFrame)
    assert df.equals(df_t)


def test_filter_smiles():
    dict_smi = dm.filter.filter_smiles(list_of_smi)

    scope = []
    catalog = []
    des = []
    for i, smi in enumerate(list_of_smi, 0):
        mol = dm.to_mol(smi)
        scope.append(",".join(dm.filter.get_props(mol, prop="Scope")))
        catalog.append(",".join(l_of_props[i]))
        des.append(",".join(l_of_descriptions[i]))
    dict_smi_t = dict(
        {
            "SMILES": list_of_smi,
            "num_of_matches": n_of_matches,
            "Catalog": catalog,
            "Description": des,
        }
    )

    assert isinstance(dict_smi, dict)
    assert isinstance(dict_smi_t, dict)
    assert dict_smi == dict_smi_t
