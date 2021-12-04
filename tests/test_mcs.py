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

    # Load/export SMARTS to check RDKit versions compatibility.
    excepted_smarts = "[#6&!R]-&!@[#6&!R]-&!@[#8&!R]-&!@[#6&R]1:&@[#6&R]:&@[#6&R]2:&@[#7&R]:&@[#6&R]:&@[#7&R]:&@[#6&R](:&@[#6&R]:&@2:&@[#6&R]:&@[#6&R]:&@1-&!@[#7&!R]-&!@[#6&!R](=&!@[#8&!R])-&!@[#6&!R]=&!@[#6&!R])-&!@[#7&!R]-&!@[#6&R]1:&@[#6&R]:&@[#6&R]:&@[#6&R]:&@[#6&R]:&@[#6&R]:&@1"
    excepted_smarts_mol = dm.from_smarts(excepted_smarts)
    excepted_smarts = dm.to_smarts(excepted_smarts_mol)

    print(smarts)

    assert smarts == excepted_smarts
