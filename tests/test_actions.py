from re import L
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


def test_add_remove_bond():
    mol = dm.to_mol("CC(=O)OC1=CC=CC=C1C(=O)O")

    tmpmol = dm.actions.remove_bond_between(mol, mol.GetAtomWithIdx(0), mol.GetAtomWithIdx(1))
    assert dm.to_inchikey(tmpmol) == "LWXRIQYZZCXYDI-UHFFFAOYSA-N"

    tmpmol = dm.actions.remove_bond_between(mol, 0, 1)
    assert dm.to_inchikey(tmpmol) == "LWXRIQYZZCXYDI-UHFFFAOYSA-N"

    tmpmol = dm.actions.remove_bond_between(mol, 0, 1, sanitize=False)
    assert dm.to_inchikey(tmpmol) == "LWXRIQYZZCXYDI-UHFFFAOYSA-N"

    tmpmol = dm.actions.remove_bond_between(
        mol, mol.GetAtomWithIdx(3), mol.GetAtomWithIdx(4), sanitize=False
    )
    tmpmol = dm.actions.add_bond_between(tmpmol, 3, 4, dm.SINGLE_BOND, sanitize=True)
    assert dm.to_inchikey(tmpmol) == "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"

    tmpmol = dm.actions.remove_bond_between(
        mol, mol.GetAtomWithIdx(3), mol.GetAtomWithIdx(4), sanitize=False
    )
    tmpmol = dm.actions.add_bond_between(
        tmpmol, tmpmol.GetAtomWithIdx(3), tmpmol.GetAtomWithIdx(4), dm.SINGLE_BOND, sanitize=True
    )
    assert dm.to_inchikey(tmpmol) == "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"


def test_update_bond():

    mol = dm.to_mol("CC(=O)OC1=CC=CC=C1C(=O)O")

    with dm.without_rdkit_log():
        assert dm.actions.update_bond(mol, 0, dm.DOUBLE_BOND) is None

        new_mol = dm.actions.update_bond(
            mol, mol.GetBondBetweenAtoms(0, 1), dm.DOUBLE_BOND, sanitize=False
        )
        assert new_mol is not None
        assert dm.to_inchikey(new_mol) == "YMTXMCUZWXKNSB-UHFFFAOYSA-N"

        new_mol = dm.actions.update_bond(mol, 0, dm.DOUBLE_BOND, sanitize=False)
        assert new_mol is not None
        assert dm.to_inchikey(new_mol) == "YMTXMCUZWXKNSB-UHFFFAOYSA-N"


def test_all_atom_join():

    mol = dm.to_mol("CC(=O)OC1=CC=CC=C1C(=O)O")
    tmpmol = dm.actions.remove_bond_between(mol, 7, 8, sanitize=False)
    mols = dm.actions.all_atom_join(tmpmol, 7, 8)
    assert dm.same_mol(mols[0], mol)

    mols = dm.actions.all_atom_join(mol, 7, 8)
    assert dm.to_inchikey(mols[0]) == "HGYLTEABMPLEAQ-UHFFFAOYSA-N"

    mols = dm.actions.all_atom_join(mol, 7, 6)
    assert dm.to_inchikey(mols[0]) == "WSZJSUJBYYBSDX-UHFFFAOYSA-N"
