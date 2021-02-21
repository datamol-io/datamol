import datamol as dm


def test_to_fp():

    smiles = "CC(=O)Oc1ccccc1C(=O)O"
    mol = dm.to_mol(smiles)

    assert dm.to_fp(mol).shape[0] == 2048
    assert dm.to_fp(mol).sum() == 29
