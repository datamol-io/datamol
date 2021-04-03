import datamol as dm


def test_freesolv():

    data = dm.data.freesolv()
    assert data.shape == (642, 4)
    assert list(data.columns) == ["iupac", "smiles", "expt", "calc"]
