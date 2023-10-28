import datamol as dm


def test_find_mcs():
    smiles_list = [
        "C=CC(=O)NCCOc1cc2ncnc(Nc3ccc(Br)cc3F)c2cc1NC(=O)C=C",
        "C=CC(=O)Nc1cc2c(Nc3ccc(F)c(Br)c3)ncnc2cc1OCCCN1CCOCC1",
        "C=CC(=O)Nc1cc2c(Nc3ccc(Br)cc3F)ncnc2cc1OCCCNC(=O)CN(C)C",
        "C=CC(=O)Nc1cc2c(Nc3ccc(Br)cc3F)ncnc2cc1OCCNC(=O)NCC",
        "C=CC(=O)Nc1cc2c(Nc3ccc(Br)cc3F)ncnc2cc1OCCNC(=O)CN(C)C",
    ]
    mols = [dm.to_mol(s) for s in smiles_list]
    smarts = dm.find_mcs(mols=mols, timeout=2)

    # NOTE(hadim): hash are different given different RDKit version
    expected_hashes = [
        # RDKit >= 2023.09
        "762f483ac10cc0f45c5aa2c790f9ef52f8dfb337",
        # RDKit <= 2023.03
        "49eff32e405d17980fad428cf4063ec52e2c5fda",
    ]

    assert dm.hash_mol(dm.from_smarts(smarts)) in expected_hashes
