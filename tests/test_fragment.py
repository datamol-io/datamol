import datamol as dm


def test_brics():
    smiles = "CCCOCc1cc(c2ncccc2)ccc1"
    mol = dm.to_mol(smiles)
    frags = dm.fragment.brics(mol)
    assert len(frags) == 9


def test_frag():
    smiles = "CCCOCc1cc(c2ncccc2)ccc1"
    mol = dm.to_mol(smiles)
    frags = dm.fragment.frag(mol)
    assert len(frags) == 9


def test_recap():
    smiles = "CCCOCc1cc(c2ncccc2)ccc1"
    mol = dm.to_mol(smiles)
    frags = dm.fragment.recap(mol)
    assert len(frags) == 3


def test_anybreak():
    smiles = "CCCOCc1cc(c2ncccc2)ccc1"
    mol = dm.to_mol(smiles)
    frags = dm.fragment.anybreak(mol)
    assert len(frags) == 9


def test_mmpa():
    smiles = "CCCOCc1cc(c2ncccc2)ccc1"
    mol = dm.to_mol(smiles)

    frags = dm.fragment.mmpa_cut(mol)
    assert len(frags) == 39
    assert "CCCOCc1cccc(-c2ccccn2)c1,C(C[*:2])[*:1],C[*:1].c1ccc(-c2cccc(CO[*:2])c2)nc1\n" in frags


def test_assemble():
    # Fragment a molecule
    smiles = "CCCOCc1cc(c2ncccc2)ccc1"
    mol = dm.to_mol(smiles)
    frags = dm.fragment.brics(mol)

    # Limit the number of fragments to work with because
    # assembling is computationally intensive.
    frags = frags[:2]

    # Assemble molecules from the list of fragments
    mols = list(dm.fragment.assemble_fragment_order(frags, max_n_mols=4))

    assert len(mols) == 4


def test_break_mol():
    smiles = "CCCOCc1cc(c2ncccc2)ccc1"
    mol = dm.to_mol(smiles)
    fragments, *_, tree = dm.fragment.break_mol(mol, randomize=False, mode="brics", returnTree=True)

    assert fragments == ["CCC", "O", "C", "c1ccncc1", "c1ccccc1"]
    assert list(tree.nodes) == [0, 1, 2, 3, 4, 5, 6, 7, 8]
    assert list(tree.edges) == [(0, 1), (0, 2), (2, 3), (2, 4), (4, 5), (4, 6), (6, 7), (6, 8)]


def test_assemble_build():
    mols = [[dm.to_mol("CCCO"), dm.to_mol("CCCCCCCO")], [dm.to_mol("CCC"), dm.to_mol("CCCCCCC")]]

    results = list(dm.fragment.build(mols))
    assert len(results) == 71

    results = list(dm.fragment.build(mols, mode="rxn"))
    assert len(results) == 0

    results = list(dm.fragment.build(mols, mode=None))
    assert len(results) == 0
