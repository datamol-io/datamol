import datamol as dm


def test_pick_atom_idx():
    smiles = "OC1=CC2CCCCC2[N:1]=C1"
    mol = dm.to_mol(smiles)

    assert isinstance(dm.actions.pick_atom_idx(mol), int)
    assert dm.actions.pick_atom_idx(mol) <= mol.GetNumAtoms()


def test_all_bond_remove():

    smiles = "OC1=CC2CCCCC2[N:1]=C1"
    mol = dm.to_mol(smiles)

    mols = dm.actions.all_bond_remove(mol)
    assert isinstance(mols, list)
