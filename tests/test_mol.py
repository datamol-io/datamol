from typing import cast

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

    mol2 = dm.atom_indices_to_mol(mol)
    for atom in mol2.GetAtoms():
        assert atom.GetIntProp("molAtomMapNumber") == atom.GetIdx()

    mol3 = dm.atom_indices_to_mol(mol, copy=True)
    for atom in mol3.GetAtoms():
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


def test_clear_mol_props():
    smiles = "Cc1ccc(CO)cc1-c1ccc2c(n1)n(C)c(=O)n2CC(C)(C)C"
    mol = dm.to_mol(smiles)

    # Set properties to the molecule
    props = dict(myname="hello", a_digit=99)
    mol = dm.set_mol_props(mol, props)

    # Check
    assert "myname" in mol.GetPropsAsDict()
    assert "a_digit" in mol.GetPropsAsDict()

    # Clear all the properties
    mol2 = dm.clear_mol_props(mol)

    # Check
    assert "myname" not in mol2.GetPropsAsDict()
    assert "a_digit" not in mol2.GetPropsAsDict()

    # Clear only a single properties
    mol3 = dm.clear_mol_props(mol, property_keys=["a_digit"])

    # Check
    assert "myname" in mol3.GetPropsAsDict()
    assert "a_digit" not in mol3.GetPropsAsDict()

    # Clear only a single properties (from a string key)
    mol3 = dm.clear_mol_props(mol, property_keys="a_digit")

    # Check
    assert "myname" in mol3.GetPropsAsDict()
    assert "a_digit" not in mol3.GetPropsAsDict()


def test_strip_mol_to_core():
    mol = dm.to_mol("CC(=O)NC1CCC2=CC(=C(C(=C2C3=CC=C(C(=O)C=C13)OC)OC)OC)OC")
    mol2 = dm.strip_mol_to_core(mol)

    assert dm.to_inchikey(mol2) == "XSMDDAFPLXKNOA-UHFFFAOYSA-N"


def test_make_scaffold_generic():
    # NOTE(hadim): pretty sure doing assert on SMARTS string is fragile and might change
    # in the future RDKit versions. So... hold and wait for it to break xD

    mol = dm.to_mol("CC(=O)NC1CCC2=CC(=C(C(=C2C3=CC=C(C(=O)C=C13)OC)OC)OC)OC")
    mol2 = dm.make_scaffold_generic(mol)
    assert (
        dm.to_smarts(mol2)
        == "[#0]-[#0](=[#0])-[#0]-[#0]1-[#0]-[#0]-[#0]2:[#0]:[#0](:[#0](:[#0](:[#0]:2-[#0]2:[#0]:[#0]:[#0](:[#0](=[#0]):[#0]:[#0]:2-1)-[#0]-[#0])-[#0]-[#0])-[#0]-[#0])-[#0]-[#0]"
    )

    mol = dm.to_mol("CC(=O)NC1CCC2=CC(=C(C(=C2C3=CC=C(C(=O)C=C13)OC)OC)OC)OC")
    mol2 = dm.make_scaffold_generic(mol, include_bonds=True)
    assert (
        dm.to_smarts(mol2)
        == "[#0][#0]([#0])[#0][#0]1[#0][#0][#0]2:[#0]:[#0](:[#0](:[#0](:[#0]:2[#0]2:[#0]:[#0]:[#0](:[#0]([#0]):[#0]:[#0]1:2)[#0][#0])[#0][#0])[#0][#0])[#0][#0]"
    )


def test_to_scaffold():
    mol = dm.to_mol("CC(=O)NC1CCC2=CC(=C(C(=C2C3=CC=C(C(=O)C=C13)OC)OC)OC)OC")
    mol2 = dm.to_scaffold_murcko(mol)
    assert dm.to_inchikey(mol2) == "XSMDDAFPLXKNOA-UHFFFAOYSA-N"

    mol = dm.to_mol("CC(=O)NC1CCC2=CC(=C(C(=C2C3=CC=C(C(=O)C=C13)OC)OC)OC)OC")
    mol2 = dm.to_scaffold_murcko(mol, make_generic=True)
    assert (
        dm.to_smarts(mol2)
        == "[#0]1-[#0]-[#0]-[#0]2:[#0]:[#0]:[#0]:[#0]:[#0]:2-[#0]2:[#0]:[#0]:[#0]:[#0](=[#0]):[#0]:[#0]:2-1"
    )


def test_compute_ring_systems():
    mol = dm.to_mol("CC(=O)NC1CCC2=CC(=C(C(=C2C3=CC=C(C(=O)C=C13)OC)OC)OC)OC")

    systems = dm.compute_ring_system(mol)
    assert len(systems) == 1
    assert len(systems[0]) == 16
    assert isinstance(systems, list)
    assert isinstance(systems[0], set)

    mol = dm.to_mol("CN1C(=O)CN=C(C2=C1C=CC(=C2)Cl)C3=CC=CC=C3")

    systems = dm.compute_ring_system(mol)
    assert len(systems) == 2
    assert len(systems[0]) == 11
    assert len(systems[1]) == 6
    assert isinstance(systems, list)
    assert isinstance(systems[0], set)


def test_hash_mol():
    # rdkit.Chem.RegistrationHash is not available for rdkit before
    # 2022.09.
    if dm.is_lower_than_current_rdkit_version("2022.09"):
        with pytest.raises(NotImplementedError):
            dm.hash_mol(dm.to_mol("CCC"))

        return

    with pytest.raises(ValueError):
        dm.hash_mol(dm.to_mol("CCC"), hash_scheme="invalid")

    data = [
        {
            "smiles": "Clc1ccc2ccccc2c1",
            "all": "334adc7dd6bf2b9384e1d5dee8a983e57bbaed80",
            "no_stereo": "6c0b4720e29e47be5062b0f13d71443a6fead14a",
            "no_tautomers": "b92bca205301bf115e5acbac38827f2fc3f3b572",
        },
        {
            "smiles": "CCCCCCCCCCCCO",
            "all": "4827ecd92b1b8b865d7f9e9c47fe8aadc5827a72",
            "no_stereo": "b2ae43fb4061ec507ffb11a61b009037500ac010",
            "no_tautomers": "31cf7325e8beb6c46dfee7ce9968252309317b1b",
        },
        {
            "smiles": "Nc1cnn(-c2ccccc2)c(=O)c1Br",
            "all": "5a2bb241fb005ba2403bbd8c25453fb714fca3b0",
            "no_stereo": "0b72a56d0b499e4e9bf58975b74ecaa86e50373a",
            "no_tautomers": "32c8840b71755d5a859963be70e8db6bd29c01e7",
        },
        {
            "smiles": "Cn1c(=O)[nH]c2[nH]c(=O)[nH]c2c1=O",
            "all": "2a1f545077bec17e3b2a6f6819dda0e55eab4442",
            "no_stereo": "8d7944f78fb61a6199f90f1505c6249b317272c0",
            "no_tautomers": "ed2f9eb1f7ce5a6ba28fad24074f975470d6937c",
        },
        {
            "smiles": "CN1CCCC(CN2c3ccccc3Sc3ccccc32)C1",
            "all": "419aef09d075056837b15da108b87f6737436de2",
            "no_stereo": "389d7084f68935d52f51fe850f72a2294e04acb3",
            "no_tautomers": "6e9142b90e22152f73c2476a5531632a32f7b8c7",
        },
        {
            "smiles": "Cc1cc(OC(=O)N(C)C)n(-c2ccccc2)n1",
            "all": "f01015c60174e7efe1fee7eb6b64804da739d104",
            "no_stereo": "274054a2de68853740504b260f6dec88cffaac28",
            "no_tautomers": "68820d4a5e08bd5476b99684192fbd73d7da64bf",
        },
        {
            "smiles": "O=[N+]([O-])c1ccc(Oc2ccc(Cl)cc2Cl)cc1",
            "all": "cfff22a4055796418862d0cdf86126118113edd0",
            "no_stereo": "f077a2f705f7e00c9206230c9c78df6641673dd8",
            "no_tautomers": "3df1b08580a5615e305cb8c5eb86f519dad6257f",
        },
        {
            "smiles": "Cc1ncc([N+](=O)[O-])n1CCO",
            "all": "973825b3a38bf44b7f264daad8c4d7631c123ec4",
            "no_stereo": "7678c0d9bb81169afd4ddf6d6f0503b126c9e1d7",
            "no_tautomers": "407c5df463de39b4605132d5f28273637952d000",
        },
        {
            "smiles": "Cc1cc(C)c2ccccc2c1",
            "all": "f2e834775b5814cf7b6d503677b5751821ea2f5c",
            "no_stereo": "1fe225caa22af11e2c32320fc6d11a8819168b92",
            "no_tautomers": "5e566a9bab54c6e733e432001428e2f0d36767cc",
        },
        {
            "smiles": "Oc1cccc2cccnc12",
            "all": "2645e016677cf9e82977b8a08a3d207e7e571f51",
            "no_stereo": "0bdcdca39192466f42f547dcd084772d13fd5f0e",
            "no_tautomers": "2f4c569ee7c73804a11758d6d16aa411e7531571",
        },
    ]

    for datum in data:
        mol = dm.to_mol(datum["smiles"])

        hash_value = dm.hash_mol(mol, hash_scheme="all")
        assert hash_value == datum["all"]

        hash_value = dm.hash_mol(mol, hash_scheme="no_stereo")
        assert hash_value == datum["no_stereo"]

        hash_value = dm.hash_mol(mol, hash_scheme="no_tautomers")
        assert hash_value == datum["no_tautomers"]


def test_to_mol_keep_hs():
    smiles = "[H]OC([H])([H])c1c([H])c([H])c(C([H])([H])[H])c(-c2nc3c(c([H])c2[H])n(C([H])([H])C(C([H])([H])[H])(C([H])([H])[H])C([H])([H])[H])c(=O)n3C([H])([H])[H])c1[H]"

    mol = dm.to_mol(smiles)
    mol = cast(dm.Mol, mol)
    assert len(mol.GetAtoms()) == 25

    mol = dm.to_mol(smiles, remove_hs=False)
    mol = cast(dm.Mol, mol)
    assert len(mol.GetAtoms()) == 50


def test_clear_atom_props():
    smiles = "Cc1ccc(CO)cc1-c1ccc2c(n1)n(C)c(=O)n2CC(C)(C)C"
    mol = dm.to_mol(smiles)

    # add the `molAtomMapNumber` property to the atoms
    mol = dm.atom_indices_to_mol(mol)

    # Check
    assert all(["molAtomMapNumber" in a.GetPropsAsDict() for a in mol.GetAtoms()])

    # Remove all the properties
    mol2 = dm.clear_atom_props(mol)

    # Check
    assert not all(["molAtomMapNumber" in a.GetPropsAsDict() for a in mol2.GetAtoms()])

    # Remove only a single property
    mol2 = dm.clear_atom_props(mol, property_keys=["molAtomMapNumber"])

    # Check
    assert not all(["molAtomMapNumber" in a.GetPropsAsDict() for a in mol2.GetAtoms()])

    # Remove only a single property from a string key
    mol3 = dm.clear_atom_props(mol, property_keys="molAtomMapNumber")

    # Check
    assert not all(["molAtomMapNumber" in a.GetPropsAsDict() for a in mol3.GetAtoms()])


def test_clear_atom_map_number():
    smiles = "Cc1ccc(CO)cc1-c1ccc2c(n1)n(C)c(=O)n2CC(C)(C)C"
    mol = dm.to_mol(smiles)

    # add the `molAtomMapNumber` property to the atoms
    mol = dm.atom_indices_to_mol(mol)

    # Check
    assert all(["molAtomMapNumber" in a.GetPropsAsDict() for a in mol.GetAtoms()])

    # Remove all the properties
    mol2 = dm.clear_atom_map_number(mol)

    # Check
    assert not all(["molAtomMapNumber" in a.GetPropsAsDict() for a in mol2.GetAtoms()])


def test_get_atom_positions():
    smiles = "[H:14][c:5]1[c:3]([c:7]([c:4]([c:6]([c:8]1[N:10]([H:18])[C:2](=[N+:11]([H:19])[H:20])[N:9]([H:16])[H:17])[H:15])[H:13])[F:1])[H:12]"
    mol = dm.to_mol(smiles, remove_hs=False)
    mol = dm.conformers.generate(mol, n_confs=1, add_hs=False)

    # Get positions
    positions_1 = dm.get_atom_positions(mol, reorder_to_atom_map_number=False)
    positions_2 = dm.get_atom_positions(mol, reorder_to_atom_map_number=True)

    # both arrays should not be equal
    assert not np.allclose(positions_1, positions_2)

    # but their sums should be
    assert np.allclose(np.sum(positions_1), np.sum(positions_2))


def test_get_atom_positions_fails():
    smiles = "CCCO"
    mol = dm.to_mol(smiles)

    with pytest.raises(ValueError):
        dm.get_atom_positions(mol)

    mol_with_conf = dm.conformers.generate(mol, n_confs=1)

    dm.get_atom_positions(mol_with_conf, reorder_to_atom_map_number=False)

    with pytest.raises(ValueError):
        dm.get_atom_positions(mol_with_conf, reorder_to_atom_map_number=True)


def test_set_atom_positions():
    smiles = "[H:14][c:5]1[c:3]([c:7]([c:4]([c:6]([c:8]1[N:10]([H:18])[C:2](=[N+:11]([H:19])[H:20])[N:9]([H:16])[H:17])[H:15])[H:13])[F:1])[H:12]"

    mol = dm.to_mol(smiles, remove_hs=False)

    positions = [
        [1.7, -6.67, 3.15],
        [0.2, 4.72, 0.78],
        [3.54, -2.64, 2.88],
        [0.43, -3.87, -0.09],
        [3.44, -0.2, 1.8],
        [0.02, -1.5, -1.0],
        [2.12, -4.54, 1.9],
        [1.5, 0.48, 0.02],
        [0.53, 7.24, 0.25],
        [1.17, 2.91, -0.85],
        [-1.22, 4.15, 2.71],
        [4.64, -3.24, 4.55],
        [-0.89, -5.43, -0.78],
        [4.52, 1.43, 2.45],
        [-1.45, -1.02, -2.48],
        [-0.15, 8.68, 1.38],
        [1.65, 7.88, -1.21],
        [2.24, 3.64, -2.15],
        [-1.96, 2.4, 3.0],
        [-2.02, 5.59, 3.71],
    ]
    positions = np.array(positions)

    # Using the atom map numbers
    mol2 = dm.set_atom_positions(
        mol=mol,
        positions=positions,
        conf_id=0,
        use_atom_map_numbers=True,
    )

    # Check the 3d flag is set
    assert mol2.GetConformers()[0].Is3D()

    # Here the ordering has been changed so only the sum will be equal
    conformer = mol2.GetConformers()[0]
    assert np.allclose(conformer.GetPositions().sum(), positions.sum())

    # Without using the atom map numbers
    # Note that in that case, the conformer will be messed up here since
    # the input positions are mapped to the atom map numbers
    mol2 = dm.set_atom_positions(
        mol=mol,
        positions=positions,
        conf_id=0,
        use_atom_map_numbers=False,
    )

    # Here the order has been kept so the positions must match
    conformer = mol2.GetConformers()[0]
    np.allclose(conformer.GetPositions(), positions)

    # Check 3d flag is not set
    positions_2d = positions
    positions_2d[:, 2] = 0

    mol3 = dm.set_atom_positions(
        mol=mol,
        positions=positions,
        conf_id=0,
        use_atom_map_numbers=True,
    )

    assert not mol3.GetConformers()[0].Is3D()


def test_set_atom_positions_fails():
    smiles = "CCCO"
    mol = dm.to_mol(smiles)

    positions = [
        [1.7, -6.67, 3.15],
        [0.2, 4.72, 0.78],
        [3.54, -2.64, 2.88],
        [0.43, -3.87, -0.09],
    ]
    positions = np.array(positions)

    # Check it works
    dm.set_atom_positions(
        mol=mol,
        positions=positions,
        conf_id=0,
        use_atom_map_numbers=False,
    )

    # Use atom map numbers but the prop is not set.
    with pytest.raises(ValueError):
        dm.set_atom_positions(
            mol=mol,
            positions=positions,
            conf_id=0,
            use_atom_map_numbers=True,
        )

    # Wrong number of dimensions of the positions array
    with pytest.raises(ValueError):
        dm.set_atom_positions(
            mol=mol,
            positions=[positions],
            conf_id=0,
            use_atom_map_numbers=True,
        )

    # Wrong shape of `positions`
    with pytest.raises(ValueError):
        dm.set_atom_positions(
            mol=mol,
            positions=positions[1:, :],
            conf_id=0,
            use_atom_map_numbers=True,
        )


def test_remove_salt():
    smiles = "CN(C)C.Cl.Cl.Br"
    mol = dm.to_mol(smiles)

    # case of success
    mol_no_salt = dm.remove_salts_solvents(mol)
    assert mol_no_salt.GetNumAtoms() == mol.GetNumAtoms() - 3

    # case to keep one salt in case the molecule is consisted by multiple salts
    smiles = "[Cl].[Ca]"
    mol = dm.to_mol(smiles)
    mol_no_salt = dm.remove_salts_solvents(mol, dont_remove_everything=True)
    assert mol_no_salt.GetNumAtoms() == 1
    mol_no_salt = dm.remove_salts_solvents(mol)
    assert mol_no_salt.GetNumAtoms() == 0

    # case salt-like atoms in the molecule are unchanged
    smiles = "CN(Br)Cl"
    mol = dm.to_mol(smiles)
    mol_no_salt = dm.remove_salts_solvents(mol)
    assert mol_no_salt.GetNumAtoms() == mol.GetNumAtoms()


def test_remove_solvent():
    smiles = "CN(C)C.CS(=O)C"
    mol = dm.to_mol(smiles)

    # case of success
    mol_no_solvent = dm.remove_salts_solvents(mol)
    assert mol_no_solvent.GetNumAtoms() == mol.GetNumAtoms() - 4

    # case solvent-like atoms in the molecule are unchanged
    smiles = "CCOc1ccccc1C(=O)O"
    mol = dm.to_mol(smiles)
    mol_no_solvent = dm.remove_salts_solvents(mol)
    assert mol_no_solvent.GetNumAtoms() == mol.GetNumAtoms()

    # case solvent is larger than molecule of interest
    smiles = (
        "CC(CCC1=CC=C(C=C1)O)NCCC2=CC(=C(C=C2)O)O.C(C1C(C(C(C(O1)OC(C(CO)O)C(C(C(=O)O)O)O)O)O)O)O"
    )
    smi_compound = "CC(CCc1ccc(O)cc1)NCCc1ccc(O)c(O)c1"
    mol = dm.to_mol(smiles)

    # largest fragment removes the wrong unit
    largest_fragment = dm.keep_largest_fragment(mol)
    assert dm.to_smiles(largest_fragment, canonical=True) != smi_compound

    # define the solvent to be removed
    mol_no_solvent = dm.remove_salts_solvents(
        mol, defn_data="C(C1C(C(C(C(O1)OC(C(CO)O)C(C(C(=O)O)O)O)O)O)O)O", defn_format="smiles"
    )
    assert dm.to_smiles(mol_no_solvent, canonical=True) == smi_compound
