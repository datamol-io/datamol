import datamol as dm


def test_freesolv():
    data = dm.data.freesolv()
    assert data.shape == (642, 4)
    assert list(data.columns) == ["iupac", "smiles", "expt", "calc"]


def test_cdk2():
    data = dm.data.cdk2()
    assert data.shape == (47, 12)
    assert list(data.columns) == [
        "smiles",
        "mol",
        "id",
        "Cluster",
        "MODEL.SOURCE",
        "MODEL.CCRATIO",
        "r_mmffld_Potential_Energy-OPLS_2005",
        "r_mmffld_RMS_Derivative-OPLS_2005",
        "b_mmffld_Minimization_Converged-OPLS_2005",
        "s_st_Chirality_1",
        "s_st_Chirality_2",
        "s_st_Chirality_3",
    ]


def test_solubility():
    data = dm.data.solubility()
    assert data.shape == (1282, 7)
    assert list(data.columns) == [
        "mol",
        "ID",
        "NAME",
        "SOL",
        "SOL_classification",
        "smiles",
        "split",
    ]


def test_chembl_drugs():
    data = dm.data.chembl_drugs()
    assert data.shape == (1935, 1)
    assert list(data.columns) == ["smiles"]


def test_chembl_samples():
    data = dm.data.chembl_samples()
    assert data.shape == (2000, 1)
    assert list(data.columns) == ["smiles"]
