import datamol as dm


def test_find_mcs_with_details():
    smiles_list = [
        "C=CC(=O)NCCOc1cc2ncnc(Nc3ccc(Br)cc3F)c2cc1NC(=O)C=C",
        "C=CC(=O)Nc1cc2c(Nc3ccc(F)c(Br)c3)ncnc2cc1OCCCN1CCOCC1",
        "C=CC(=O)Nc1cc2c(Nc3ccc(Br)cc3F)ncnc2cc1OCCCNC(=O)CN(C)C",
        "C=CC(=O)Nc1cc2c(Nc3ccc(Br)cc3F)ncnc2cc1OCCNC(=O)NCC",
        "C=CC(=O)Nc1cc2c(Nc3ccc(Br)cc3F)ncnc2cc1OCCNC(=O)CN(C)C",
    ]
    mols = [dm.to_mol(s) for s in smiles_list]
    mcs = dm.find_mcs_with_details(mols=mols, timeout=1)

    excepted_smarts = "[#6]-[#6]-[#8]-[#6]1:[#6]:[#6]2:[#7]:[#6]:[#7]:[#6](:[#6]:2:[#6]:[#6]:1-[#7]-[#6](=[#8])-[#6]=[#6])-[#7]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1"

    assert mcs.smartsString == excepted_smarts


def test_find_mcs():
    smiles_list = [
        "C=CC(=O)NCCOc1cc2ncnc(Nc3ccc(Br)cc3F)c2cc1NC(=O)C=C",
        "C=CC(=O)Nc1cc2c(Nc3ccc(F)c(Br)c3)ncnc2cc1OCCCN1CCOCC1",
        "C=CC(=O)Nc1cc2c(Nc3ccc(Br)cc3F)ncnc2cc1OCCCNC(=O)CN(C)C",
        "C=CC(=O)Nc1cc2c(Nc3ccc(Br)cc3F)ncnc2cc1OCCNC(=O)NCC",
        "C=CC(=O)Nc1cc2c(Nc3ccc(Br)cc3F)ncnc2cc1OCCNC(=O)CN(C)C",
    ]
    mols = [dm.to_mol(s) for s in smiles_list]
    smarts = dm.find_mcs(mols=mols, timeout=1)

    excepted_smarts = "[#6]-[#6]-[#8]-[#6]1:[#6]:[#6]2:[#7]:[#6]:[#7]:[#6](:[#6]:2:[#6]:[#6]:1-[#7]-[#6](=[#8])-[#6]=[#6])-[#7]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1"

    assert smarts == excepted_smarts
