import datamol as dm


def test_without_rdkit_log(capfd):
    """NOTE(hadim): unittest.TestCase does not workk with `unittest.TestCase`."""
    smiles = "fake_smiles"
    mol = dm.to_mol(smiles)
    _, err = capfd.readouterr()
    assert "SMILES Parse Error: syntax error while parsing: fake_smiles" in err

    with dm.without_rdkit_log():
        smiles = "fake_smiles"
        mol = dm.to_mol(smiles)
        _, err = capfd.readouterr()
        assert err == ""
