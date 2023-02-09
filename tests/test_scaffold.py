import datamol as dm


def test_fuzzy_scaffolding():
    smiles = [
        "Cc1ccc(NC(=O)Cn2cccn2)c(Br)c1",
        "COc1ccc(OC(C)C(=O)N=c2sccn2C)cc1",
        "CC(NC(=O)CSCc1cccs1)C1CCCO1",
        "CC1CCCCN1C(=O)CN1CCC[C@@H](N)C1",
        "CCC(CC)COC(=O)[C@H](C)N[P@](=O)(OC[C@H]1O[C@](C#N)([C@H](O)[C@@H]1O)C1=CC=C2N1N=CN=C2N)OC1=CC=CC=C1",  # no way this one (Remdesivir) is in the db
    ]

    mols = [dm.to_mol(s) for s in smiles]
    all_scaffolds, scf2infos, scf2groups = dm.scaffold.fuzzy_scaffolding(mols)

    assert scf2infos.keys() == scf2groups.keys()
    assert len(all_scaffolds) == 5

    # NOTE(hadim): different version of rdkit (2020.09 vs 2021.03) returns
    # different SMILES here.
    # assert "O=C(CN1CCC[C@@H]([*:1])C1)N1CCCCC1[*:2]" in all_scaffolds
    # assert "O=C(CSCc1cccs1)NC(C1CCCO1)[*:1]" in all_scaffolds
    # assert "O=C(N=c1sccn1[*:1])C(Oc1ccc([*:3])cc1)[*:2]" in all_scaffolds
