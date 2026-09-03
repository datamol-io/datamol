import datamol as dm


def test_fuzzy_scaffolding():
    smiles = [
        "Cc1ccc(NC(=O)Cn2cccn2)c(Br)c1",
        "COc1ccc(OC(C)C(=O)N=c2sccn2C)cc1",
        "CC(NC(=O)CSCc1cccs1)C1CCCO1",
        "CC1CCCCN1C(=O)CN1CCC[C@@H](N)C1",
        "CCC(CC)COC(=O)[C@H](C)N[P@](=O)(OC[C@H]1O[C@](C#N)([C@H](O)[C@@H]1O)C1=CC=C2N1N=CN=C2N)OC1=CC=CC=C1",  # no way this one (Remdesivir) is in the db
        "COc1ccc(OC(C)C(=O)N=c2sccn2C)cc1",
    ]

    mols = [dm.to_mol(s) for s in smiles]

    # NOTE(hadim): different version of rdkit (2020.09 vs 2021.03) returns
    # different SMILES here.
    # assert "O=C(CN1CCC[C@@H]([*:1])C1)N1CCCCC1[*:2]" in all_scaffolds
    # assert "O=C(CSCc1cccs1)NC(C1CCCO1)[*:1]" in all_scaffolds
    # assert "O=C(N=c1sccn1[*:1])C(Oc1ccc([*:3])cc1)[*:2]" in all_scaffolds

    all_scaffolds, df_scf2infos, df_scf2groups = dm.scaffold.fuzzy_scaffolding(mols)

    assert len(all_scaffolds) == 5
    assert len(df_scf2infos.columns) == 3

    # because we are returning the output for each scf
    # these should be the same
    assert len(df_scf2infos.index) == len(df_scf2groups.index)
    assert list(df_scf2infos["scf"]) == list(df_scf2groups["scf"])

    # mere coincidence that scf2infos and scf2groups for the columns have the
    # the same length. the reason there are 3 not two is because it could have
    # extra columns where a cell may have none values.
    assert len(df_scf2groups.columns) == 3


def test_fuzzy_scaffolding_with_enforce_subs():
    """Test that enforce_subs accepts string (SMILES/SMARTS) arguments.

    Regression test for https://github.com/datamol-io/datamol/issues/119
    where passing strings to enforce_subs caused a TypeError in
    GetSubstructMatch which expects an RDKit Mol object.
    """
    smiles = [
        "Cc1ccc(NC(=O)Cn2cccn2)c(Br)c1",
        "COc1ccc(OC(C)C(=O)N=c2sccn2C)cc1",
        "CC(NC(=O)CSCc1cccs1)C1CCCO1",
        "CC1CCCCN1C(=O)CN1CCC[C@@H](N)C1",
        "COc1ccc(OC(C)C(=O)N=c2sccn2C)cc1",
    ]

    mols = [dm.to_mol(s) for s in smiles]

    # Pass enforce_subs as strings (SMILES) — should not raise
    all_scaffolds, df_scf2infos, df_scf2groups = dm.scaffold.fuzzy_scaffolding(
        mols, enforce_subs=["c1ccccc1", "C1CCCO1"]
    )

    assert isinstance(all_scaffolds, set)
    assert len(all_scaffolds) >= 1
    assert len(df_scf2infos.index) == len(df_scf2groups.index)

    # Pass enforce_subs as SMARTS strings — should also not raise
    all_scaffolds2, _, _ = dm.scaffold.fuzzy_scaffolding(
        mols, enforce_subs=["[cR1]1[cR1][cR1][cR1][cR1][cR1]1"]
    )

    assert isinstance(all_scaffolds2, set)
