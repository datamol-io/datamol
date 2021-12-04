import datamol as dm


def test_match_substructure():
    mol1 = dm.to_mol("CC(=O)OC1=CC=CC=C1C(=O)O")
    mol2 = dm.to_mol("CCN(CC)CC(=O)CC(C)NC1=C2C=CC(=CC2=NC=C1)Cl")

    query1 = dm.from_smarts("[C;H0](=O)")
    query2 = dm.to_mol("CN(C)")

    # Test multiple scenarios

    dm.viz.match_substructure(
        mols=[mol1, mol2],
        queries=[query1, query2],
        highlight_bonds=True,
        use_svg=True,
    )
    dm.viz.match_substructure(
        mols=mol1,
        queries=[query1, query2],
        highlight_bonds=True,
        use_svg=True,
    )
    dm.viz.match_substructure(
        mols=[mol1, mol2],
        queries=query1,
        highlight_bonds=False,
        use_svg=False,
    )
    dm.viz.match_substructure(
        mols=mol1,
        queries=query2,
        highlight_bonds=True,
        use_svg=False,
    )
