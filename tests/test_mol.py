import pytest

import itertools

from rdkit import Chem

import datamol as dm
import numpy as np


def test_to_mol():
    smiles = "O=C(C)Oc1ccccc1C(=O)O"
    mol = dm.to_mol(smiles)
    assert mol.GetNumAtoms() == 13

    smiles = "O=C(C)Oc1ccccc1C(=O)O"
    mol = dm.to_mol(smiles, add_hs=True)
    assert mol.GetNumAtoms() == 21

    smiles = "fake_smiles"
    mol = dm.to_mol(smiles)
    assert mol is None


def test_reorder_atoms():
    smiles = "c1ccc(C(=O)O)c(c1)OC(=O)C"
    mol = dm.to_mol(smiles, add_hs=False, explicit_only=False)

    orders = [a.GetAtomicNum() for a in mol.GetAtoms()]
    assert orders == [6, 6, 6, 6, 6, 8, 8, 6, 6, 8, 6, 8, 6]

    mol = dm.reorder_atoms(mol)
    orders = [a.GetAtomicNum() for a in mol.GetAtoms()]
    assert orders == [6, 8, 8, 8, 6, 6, 6, 6, 8, 6, 6, 6, 6]


def test_randomize_atoms():
    smiles = "c1ccc(C(=O)O)c(c1)OC(=O)C"
    mol = dm.to_mol(smiles)
    orders = [a.GetAtomicNum() for a in mol.GetAtoms()]

    randomized_mol = dm.randomize_atoms(mol)
    randomized_orders = [a.GetAtomicNum() for a in randomized_mol.GetAtoms()]

    assert sum(orders) == sum(randomized_orders)


def test_to_neutral():

    smiles = "[NH4+]"
    mol = dm.to_mol(smiles, add_hs=False, explicit_only=False)

    smiles = dm.to_smiles(dm.to_neutral(mol))
    assert smiles == "[NH4]"

    smiles = "O=C(c1ccccc1)[O-]"
    mol = dm.to_mol(smiles, add_hs=False, explicit_only=False)
    uncharged_mol = dm.to_neutral(mol)
    assert sum([a.GetFormalCharge() for a in uncharged_mol.GetAtoms()]) == 0


def test_sanitize():
    smiles = "CC(=O)Oc1ccccc1C(=O)O"
    mol = dm.to_mol(smiles, sanitize=False)
    mol = dm.sanitize_mol(mol, charge_neutral=True)
    assert dm.to_smiles(mol) == "CC(=O)Oc1ccccc1C(=O)O"

    mol = dm.sanitize_mol(None, charge_neutral=True)
    assert mol is None

    smiles_list = (
        "CC.[H][N:1]1(C)=CC(O)=CC2CCCCC12",  # broken
        "O=c1ccc2ccccc2n1",  # sanitize
        "Cc1nnnn1C",  # none
        "CCc1ccc2nc(=O)c(cc2c1)Cc1nnnn1C1CCCCC1",  # sanitize
        "c1cnc2cc3ccnc3cc12",  # none
        "c1cc2cc3ccnc3cc2n1",  # none
        "O=c1ccnc(c1)-c1cnc2cc3ccnc3cc12",  # sanitize
        "O=c1ccnc(c1)-c1cc1",  # broken
    )

    # check sanitize_mol
    assert dm.to_mol(smiles_list[1]) is None
    assert dm.to_mol(smiles_list[2]) is not None
    assert dm.sanitize_mol(None) is None
    assert dm.sanitize_mol(dm.to_mol(smiles_list[0], sanitize=False)) is None
    assert dm.sanitize_mol(dm.to_mol(smiles_list[1], sanitize=False)) is not None
    assert dm.sanitize_mol(dm.to_mol(smiles_list[2], sanitize=False)) is not None

    mol_2 = dm.sanitize_mol(dm.to_mol(smiles_list[1], sanitize=False))
    assert dm.to_smiles(mol_2) == dm.sanitize_smiles("O=c1ccc2ccccc2[nH]1")

    fixed_smiles = [dm.sanitize_smiles(smiles) for smiles in smiles_list]
    assert len([x for x in fixed_smiles if x is not None]) == 6

    smiles = "CCCN(=O)=O"
    mol = dm.to_mol(smiles, sanitize=False)
    sane_mol = dm.sanitize_mol(mol)
    assert dm.to_smiles(sane_mol) == "CCC[N+](=O)[O-]"


def test_sanitize_first():

    smiles = ["fake_smiles", "CC(=O)Oc1ccccc1C(=O)O"]
    mols = [dm.to_mol(s) for s in smiles]
    mol = dm.sanitize_first(mols)
    assert dm.to_smiles(mol) == "CC(=O)Oc1ccccc1C(=O)O"


def test_standardize_mol():
    sm = "[Na]OC1=CC2CCCCC2N=C1"
    sm_standard = dm.to_smiles(dm.standardize_smiles(sm))
    standard_mol = dm.standardize_mol(dm.to_mol(sm), disconnect_metals=True, uncharge=True)
    mol_standard = dm.to_smiles(Chem.MolToSmiles(standard_mol))
    assert sm_standard == mol_standard


def test_fix_valence():
    sm = "Cl.[H][N:1]1=CC(O)=CC2CCCCC12"
    mol = Chem.MolFromSmiles(sm, sanitize=False)
    mol.UpdatePropertyCache(False)
    mol_copy = dm.copy_mol(mol)

    nitrogen_atom = [a for a in mol.GetAtoms() if a.GetAtomMapNum() == 1][0]
    nitrogen_valence = nitrogen_atom.GetExplicitValence()
    assert dm.incorrect_valence(nitrogen_atom, True)

    fixed_mol = dm.fix_valence_charge(mol, inplace=False)
    assert dm.to_mol(Chem.MolToSmiles(fixed_mol)) is not None

    # expect nitrogen atom to still be incorrect
    assert dm.incorrect_valence(nitrogen_atom, True)

    # in place fix
    fixed_mol = dm.fix_valence_charge(mol, inplace=True)
    # nitrogen should be charged positively if this was fixed.
    assert nitrogen_atom.GetFormalCharge() == 1

    fixed_mol2 = dm.fix_valence(mol_copy)
    fixed_nitrogen_atom = [a for a in fixed_mol2.GetAtoms() if a.GetAtomMapNum() == 1][0]
    assert fixed_nitrogen_atom.GetExplicitValence() < nitrogen_valence

    # mol should be fixed
    assert dm.to_mol(Chem.MolToSmiles(fixed_mol2)) is not None


def test_adjust_singleton():
    sm = "Cl.[N:1]1=CC(O)=CC2CCCCC12.CC.C"
    mol = dm.to_mol(sm)
    fixed_mol = dm.adjust_singleton(mol)
    assert len(Chem.rdmolops.GetMolFrags(fixed_mol)) == 2
    assert fixed_mol.HasSubstructMatch(Chem.MolFromSmiles("CC"))  # assert ethyl is there


def test_fixmol():
    sm = "C.Cl.CC.[H][N:1]1(C)=CC(O)=CC2CCCCC12"
    mol = Chem.MolFromSmiles(sm, sanitize=False)
    # mol.UpdatePropertyCache(False)
    # Chem.Kekulize(mol)
    res = dm.fix_mol(mol, n_iter=1)  # copy by default

    # should still be invalid in term of valence for nitrogen
    assert not dm.incorrect_valence(res)

    res2 = dm.fix_mol(mol, n_iter=2)
    # not expecting difference between res2 and res3
    assert Chem.MolToSmiles(res) == Chem.MolToSmiles(res2)

    # only largest expected_here
    res_largest = dm.fix_mol(mol, largest_only=True)

    dm.fix_mol(mol, remove_singleton=True, largest_only=True)
    assert len(Chem.rdmolops.GetMolFrags(res_largest)) == 1

    expected_largest_fix = dm.standardize_smiles("OC1=CC2CCCCC2[N:1]=C1")
    assert dm.standardize_smiles(Chem.MolToSmiles(res_largest)) == expected_largest_fix

    res_no_singleton = dm.fix_mol(mol, n_iter=2, remove_singleton=True)
    assert len(Chem.rdmolops.GetMolFrags(res_largest)) == 1
    assert len(Chem.rdmolops.GetMolFrags(res_no_singleton)) == 2


def test_dative_bond():
    smis = "CC1=CC=CC(=C1N\\2O[Co]3(ON(\\C=[N]3\\C4=C(C)C=CC=C4C)C5=C(C)C=CC=C5C)[N](=C2)\\C6=C(C)C=CC=C6C)C"
    expected_result = (
        "CC1=CC=CC(C)=C1N1C=N(C2=C(C)C=CC=C2C)->[Co]2(<-N(C3=C(C)C=CC=C3C)=CN(C3=C(C)C=CC=C3C)O2)O1"
    )

    assert dm.is_transition_metal(Chem.Atom("Co"))

    # sodium is not a transition metal
    assert not dm.is_transition_metal(Chem.Atom("Na"))

    mol = dm.set_dative_bonds(Chem.MolFromSmiles(smis, sanitize=False))
    assert Chem.MolToSmiles(mol) == expected_result
    assert dm.to_mol(Chem.MolToSmiles(mol)) is not None


def test_copy_mol():
    mol = dm.to_mol("OC1=CC2CCCCC2[N:1]=C1")
    new_mol = dm.copy_mol(mol)

    assert dm.to_smiles(mol) == dm.to_smiles(new_mol)


def test_set_mol_props():
    mol = dm.to_mol("CCC")

    props = {}
    props["number"] = 55
    props["float"] = 5.555
    props["string"] = "hello"
    props["something_else"] = type(int)

    dm.set_mol_props(mol, props)

    mol_props = mol.GetPropsAsDict()
    assert mol_props["number"] == props["number"]
    assert mol_props["float"] == props["float"]
    assert mol_props["string"] == props["string"]
    assert mol_props["something_else"] == str(props["something_else"])

    dm.set_mol_props(mol, props, copy=True)


def test_set_mol_props_overflow():

    mol = dm.to_mol("CCC")
    dm.set_mol_props(mol, dict(hello=661440088496))
    assert mol.GetPropsAsDict() == {"hello": 661440088496}


def test_copy_mol_props():
    source = dm.to_mol("CCC")
    destination = dm.to_mol("CC")

    props = {}
    props["bool"] = True
    props["number"] = 55
    props["float"] = 5.555
    props["string"] = "hello"
    props["something_else"] = type(int)

    dm.set_mol_props(source, props)

    dm.copy_mol_props(source, destination)

    assert destination.GetPropsAsDict() == source.GetPropsAsDict()


def test_atom_indices_to_mol():
    mol: dm.Mol = dm.to_mol("OC1=CC2CCCCC2[N:1]=C1")

    dm.atom_indices_to_mol(mol)
    for atom in mol.GetAtoms():
        assert atom.GetIntProp("molAtomMapNumber") == atom.GetIdx()

    dm.atom_indices_to_mol(mol, copy=True)
    for atom in mol.GetAtoms():
        assert atom.GetIntProp("molAtomMapNumber") == atom.GetIdx()


def test_to_mol_bad_input():
    with pytest.raises(ValueError):
        dm.to_mol(555)


def test_to_mol_input_mol():
    mol = dm.to_mol(dm.to_mol("CCC"))
    assert dm.to_smiles(mol) == "CCC"


def test_to_mol_ordered():
    s1 = dm.to_smiles(dm.to_mol("C1(O)=CC2CCCCC2[N:1]=C1", ordered=False), canonical=False)
    s2 = dm.to_smiles(dm.to_mol("C1(O)=CC2CCCCC2[N:1]=C1", ordered=True), canonical=False)
    assert s1 != s2


def test_to_mol_kekulize():
    smiles = "C1=CC=CN=C1"
    dm.to_mol(smiles, kekulize=True)


def test_reorder_atoms_without_atoms():
    mol1 = dm.to_mol("")
    mol2 = dm.reorder_atoms(mol1)
    assert dm.to_smiles(mol1) == dm.to_smiles(mol2)


def test_randomize_atoms_without_atoms():
    mol1 = dm.to_mol("")
    mol2 = dm.randomize_atoms(mol1)
    assert dm.to_smiles(mol1) == dm.to_smiles(mol2)


def test_to_neutral_without_atoms():
    assert dm.to_neutral(None) is None


def test_sanitize_smiles_none():
    assert dm.sanitize_smiles(444) is None


def test_standardize_smiles_tautomer():
    smiles = "C1=CC=CN=C1"
    std_smiles = dm.standardize_smiles(smiles, tautomer=True)
    assert "c1ccncc1" == std_smiles


def test_decrease_bond():

    smiles = "C=CCC#C"
    mol = dm.to_mol(smiles)

    single_bond = mol.GetBondBetweenAtoms(1, 2)
    double_bond = mol.GetBondBetweenAtoms(0, 1)
    triple_bond = mol.GetBondBetweenAtoms(3, 4)

    assert dm.decrease_bond(single_bond) is None
    assert dm.decrease_bond(double_bond) == dm.SINGLE_BOND
    assert dm.decrease_bond(triple_bond) == dm.DOUBLE_BOND


def test_replace_dummies_atoms():
    smiles = "C=CCC*"
    mol = dm.to_mol(smiles)

    mol2 = dm.replace_dummies_atoms(mol)
    assert dm.to_smiles(mol2) == "C=CCCC"


def test_keep_largest_frag():
    mol = dm.to_mol("Cl.[N:1]1=CC(O)=CC2CCCCC12.CC.C")
    frag = dm.keep_largest_fragment(mol)
    assert dm.to_smiles(frag) == "OC1=CC2CCCCC2[N:1]=C1"


def test_sanitize_mol_keep_props_and_conformers():
    # Generate a mol with props and a conformer
    props = dict(test_int=588, test_str="hello")
    smiles = "CCC[N+](=O)[O-]"

    mol = dm.to_mol(smiles)
    mol = dm.set_mol_props(mol, props)
    mol = dm.conformers.generate(mol, n_confs=1)
    pos = mol.GetConformer().GetPositions()

    # Sanitize
    sane_mol = dm.sanitize_mol(mol)

    # Check properties
    assert sane_mol.GetPropsAsDict() == props

    # Check conformer
    conf = sane_mol.GetConformer()
    assert sane_mol.GetNumConformers() == 1
    assert conf.Is3D()
    np.testing.assert_almost_equal(conf.GetPositions(), pos, decimal=4)


def test_sanitize_mol_multiple_conformers_warning(caplog):
    # Generate a mol with props and a conformer
    smiles = "CCC[N+](=O)[O-]"

    mol = dm.to_mol(smiles)
    mol = dm.conformers.generate(mol, n_confs=10)

    # Check warning log
    dm.sanitize_mol(mol)
    assert "WARNING" != caplog.text
    assert (
        "The molecule contains multiple conformers. Only the first one will be preserved."
        in caplog.text
    )


def test_sanitize_mol_multiple_conformers_no_warning(caplog):

    # Generate a mol with props and a conformer
    smiles = "CCC[N+](=O)[O-]"

    mol = dm.to_mol(smiles)
    mol = dm.conformers.generate(mol, n_confs=10)

    # Check no warning log
    dm.sanitize_mol(mol, verbose=False)
    assert caplog.text == ""


def test_same_mol():
    mol1 = dm.to_mol("CC(=O)Oc1ccccc1C(=O)O")
    mol2 = dm.to_mol("C1OC1CC")

    assert dm.same_mol(mol1, mol2) is False
    assert dm.same_mol(mol1, mol1) is True
    assert dm.same_mol(mol2, mol2) is True
    assert dm.same_mol(None, mol2) is False
    assert dm.same_mol(mol1, None) is False
    assert dm.same_mol(None, None) is False

    # check with not sane molecule
    mol1 = dm.to_mol("c1ccccc1")
    mol2 = dm.to_mol("C1=CC=CC=C1", sanitize=False)  # not sane benzene but still a benzene molecule

    assert dm.same_mol(mol1, mol2) is True
    assert dm.same_mol(mol1, mol1) is True
    assert dm.same_mol(mol2, mol2) is True
    assert dm.same_mol(None, mol2) is False
    assert dm.same_mol(mol1, None) is False
    assert dm.same_mol(None, None) is False

    mol1 = dm.to_mol("N=C(N)O")
    mol2 = dm.to_mol("NC(N)=O")

    assert dm.same_mol(mol1, mol2) is True
    assert dm.same_mol(mol1, mol2, use_non_standard_inchikey=True) is False


def test_atom_list_to_bond():
    mol = dm.to_mol("CC(=O)OC1=CC=CC=C1C(=O)O")

    bond_indices = dm.atom_list_to_bond(mol, atom_indices=[3, 4, 5, 6, 10, 11], bond_as_idx=True)
    assert bond_indices == [3, 4, 5, 10]

    bonds = dm.atom_list_to_bond(mol, atom_indices=[3, 4, 5, 6, 10, 11], bond_as_idx=False)
    assert [b.GetIdx() for b in bonds] == [3, 4, 5, 10]


def test_protect_atoms():
    mol = dm.to_mol("CC(=O)Oc1cnccc1C(O)=O")

    # carbonyl group
    query = dm.from_smarts("[$([CX3]=[OX1]),$([CX3+]-[OX1-])]")

    # find id of single nitrogen
    atom_ids = [at.GetIdx() for at in mol.GetAtoms() if at.GetSymbol() == "N"]  # type: ignore
    protected_atoms = itertools.chain(*mol.GetSubstructMatches(query))  # type: ignore
    protected_atoms = set(list(protected_atoms) + atom_ids)
    new_mol = dm.protect_atoms(mol, substruct=query, atoms=atom_ids)
    computed_protect_atoms = set(
        [
            at.GetIdx()
            for at in new_mol.GetAtoms()
            if int(at.GetPropsAsDict().get("_protected", 0)) == 1
        ]
    )
    assert protected_atoms == computed_protect_atoms


def test_add_remove_hs():

    smiles = "OC1=CC2CCCCC2[N:1]=C1"

    mol = dm.to_mol(smiles)

    mol2 = dm.add_hs(mol)
    assert (
        dm.to_smiles(mol2)
        == "[H]OC1=C([H])C2([H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])C2([H])[N:1]=C1[H]"
    )

    mol3 = dm.remove_hs(mol)
    assert dm.to_smiles(mol3) == smiles


def test_unique_id():
    smiles_list = [
        # mols coming from FreeSolv
        "c1ccc2c(c1)Cc3ccccc3C2",
        "CCCC/C=C/C",
        "Cc1c[nH]cn1",
        "Cc1ccc(nc1)C",
        "Cc1cccs1",
        "C(CO[N+](=O)[O-])CO[N+](=O)[O-]",
        "C(=C(F)F)(C(F)(F)F)F",
        "CCOP(=S)(OCC)SCSP(=S)(OCC)OCC",
        "CC1=CC(=O)[C@@H](CC1)C(C)C",
        "Cc1cc(ccc1Cl)O",
        "CC(C)O",
        "C",
        "CCCC(=O)C",
        "CCCCCCCC=O",
        "CCc1cccc2c1cccc2",
        # Two cherry picked tautomeric forms
        # of the same molecule.
        "N=C(N)O",
        "NC(N)=O",
        # an invalid molecule
        "C(C)(C)(C)(C)(C)",
        # a molecule with a dummy atom
        "CCCC*",
    ]

    mols = [dm.to_mol(s) for s in smiles_list]
    ids = [dm.unique_id(m) for m in mols]
    assert ids == [
        "2f0a84f817a369e79a2b44085ebe3a27",
        "3dc94b495eba316c2df4563457c77842",
        "df2c2dcbe7262f02f371d5523248e3d6",
        "1d0bc775e2b6ab0a077548a26ce9586f",
        "6cb8a819d03a38a6ed4771d1bbfd899c",
        "ae2683fabd3b5c31ac0554c213e1f37a",
        "c77fd251ff9b35379e8b08f1cf870391",
        "1a239153353f30ec3067af0461a5bff3",
        "f5bcdce1e7726f01abdc8f30f627ebcf",
        "4121e41877275378c2ae501815ab2876",
        "a6e46ea133142499b60f73a6d3e94248",
        "cbaaa6452f9cc0649cc9920a698c6f5d",
        "70572a2d32cf4707138a6ef81cfd7253",
        "180d723f0fbcc955a975d740fbc91adf",
        "11e555af02fa8342860752227c0ddc9e",
        # The two IDS of the tautomers are indeed different.
        "3a04e3e47b86287b035b4affe84dff91",
        "943d4548ef05069b79e929ed140ce905",
        # unsupported syntax should return None
        None,
        None,
    ]

    assert dm.unique_id(None) is None
