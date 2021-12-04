import datamol as dm


def test_align_2d_coordinates():
    smiles_list = [
        "COCCOc1cc2ncnc(Nc3ccc(Br)cc3F)c2cc1NC(=O)/C=C/CN(C)C",
        "C=CC(=O)Nc1cc2c(Nc3ccc(Br)cc3F)ncnc2cc1OCCNCC(C)=O",
        "C=CC(=O)Nc1cc2c(Nc3cccc(Br)c3)ncnc2cc1OCCCN1CCOCC1",
        "C=CC(=O)Nc1cc2c(Nc3cc(Cl)c(Cl)cc3F)ncnc2cc1OCCNCC(C)=O",
        "C=CC(=O)Nc1cc2c(Nc3ccc(Br)cc3F)ncnc2cc1OCCCNC(C)=O",
        "C=CC(=O)Nc1cc2c(Nc3ccc(Br)cc3F)ncnc2cc1OCCN",
        "C=CC(=O)Nc1cc2c(Nc3ccc(Br)cc3F)ncnc2cc1OCCNC(N)=O",
        "C=CC(=O)Nc1cc2c(Nc3cc(Cl)c(Cl)cc3F)ncnc2cc1OCCNC(C)=O",
        "C=CC(=O)Nc1cc2c(Nc3ccc(Br)cc3F)ncnc2cc1OCCNC(=O)C(F)F",
    ]
    mols = [dm.to_mol(s) for s in smiles_list]

    dm.viz.utils.align_2d_coordinates(mols)


def test_align_2d_coordinates_with_pattern():

    smiles_list = [
        "COCCOc1cc2ncnc(Nc3ccc(Br)cc3F)c2cc1NC(=O)/C=C/CN(C)C",
        "C=CC(=O)Nc1cc2c(Nc3ccc(Br)cc3F)ncnc2cc1OCCNCC(C)=O",
        "C=CC(=O)Nc1cc2c(Nc3cccc(Br)c3)ncnc2cc1OCCCN1CCOCC1",
        "C=CC(=O)Nc1cc2c(Nc3cc(Cl)c(Cl)cc3F)ncnc2cc1OCCNCC(C)=O",
        "C=CC(=O)Nc1cc2c(Nc3ccc(Br)cc3F)ncnc2cc1OCCCNC(C)=O",
        "C=CC(=O)Nc1cc2c(Nc3ccc(Br)cc3F)ncnc2cc1OCCN",
        "C=CC(=O)Nc1cc2c(Nc3ccc(Br)cc3F)ncnc2cc1OCCNC(N)=O",
        "C=CC(=O)Nc1cc2c(Nc3cc(Cl)c(Cl)cc3F)ncnc2cc1OCCNC(C)=O",
        "C=CC(=O)Nc1cc2c(Nc3ccc(Br)cc3F)ncnc2cc1OCCNC(=O)C(F)F",
    ]
    mols = [dm.to_mol(s) for s in smiles_list]

    pattern = dm.from_smarts(
        "[#6]=[#6]-[#6](=[#8])-[#7]-[#6]1:[#6]:[#6]2:[#6](-[#7]-[#6]3:[#6]:[#6]:[#6]:[#6]:[#6]:3):[#7]:[#6]:[#7]:[#6]:2:[#6]:[#6]:1-[#8]-[#6]-[#6]"
    )

    dm.viz.utils.align_2d_coordinates(mols, pattern=pattern)
