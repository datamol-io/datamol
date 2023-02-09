import datamol as dm
import unittest as ut


def test_to_mol():
    smiles = "O=C(C)Oc1ccccc1C(=O)O"
    mol = dm.to_mol(smiles)
    graph = dm.to_graph(mol)

    assert list(graph.nodes) == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    assert list(graph.edges) == [
        (0, 1),
        (1, 2),
        (1, 3),
        (3, 4),
        (4, 5),
        (4, 9),
        (5, 6),
        (6, 7),
        (7, 8),
        (8, 9),
        (9, 10),
        (10, 11),
        (10, 12),
    ]


def test_get_all_path_between():
    smiles = "c1cc2cccccc2c1"
    mol = dm.to_mol(smiles)

    all_paths = dm.get_all_path_between(mol, 8, 4, ignore_cycle_basis=False)
    assert all_paths == [[8, 2, 3, 4], [8, 7, 6, 5, 4], [8, 9, 0, 1, 2, 3, 4]]

    all_paths = dm.get_all_path_between(mol, 8, 4, ignore_cycle_basis=True)
    assert all_paths == [[8, 2, 3, 4], [8, 7, 6, 5, 4]]


class Test_match_molecular_graphs(ut.TestCase):
    def test_ring(self):
        mol1 = dm.to_mol("C1CCCCC1", ordered=False)
        mol2 = dm.to_mol("C1=CC=CC=C1", ordered=True)
        mol3 = dm.to_mol("C1CCOCC1", ordered=False)
        mol4 = dm.to_mol("C1=CCOCC1", ordered=True)
        mol5 = dm.to_mol("C1=CCCCC1", ordered=True)
        mol6 = dm.to_mol("C1CCCCC1", ordered=True, add_hs=True)

        # 12 matches on rings (6 rotations * 2 reflections)
        matches = dm.match_molecular_graphs(
            mol1, mol1, match_atoms_on="atomic_num", match_bonds_on=["bond_type"]
        )
        self.assertEqual(len(matches), 12)

        # 12 matches on rings (6 rotations * 2 reflections)
        matches = dm.match_molecular_graphs(
            mol1, mol2, match_atoms_on="atomic_num", match_bonds_on=[]
        )
        self.assertEqual(len(matches), 12)

        # 12 matches on rings (6 rotations * 2 reflections)
        matches = dm.match_molecular_graphs(
            mol2, mol2, match_atoms_on=[], match_bonds_on=["bond_type"]
        )
        self.assertEqual(len(matches), 12)

        # 0 matches on rings of different bond types
        matches = dm.match_molecular_graphs(
            mol1, mol2, match_atoms_on=[], match_bonds_on=["bond_type"]
        )
        self.assertEqual(len(matches), 0)

        # Matches with Oxygen atom
        matches = dm.match_molecular_graphs(
            mol1, mol3, match_atoms_on=[], match_bonds_on=["bond_type"]
        )
        self.assertEqual(len(matches), 12)

        # No matching due to atom type
        matches = dm.match_molecular_graphs(
            mol1, mol3, match_atoms_on=["atomic_num"], match_bonds_on=[]
        )
        self.assertEqual(len(matches), 0)

        # 2 reflections can match
        matches = dm.match_molecular_graphs(
            mol3, mol4, match_atoms_on=["atomic_num"], match_bonds_on=[]
        )
        self.assertEqual(len(matches), 2)

        # No matching due to different bonds
        matches = dm.match_molecular_graphs(
            mol3, mol4, match_atoms_on=["atomic_num"], match_bonds_on=["bond_type"]
        )
        self.assertEqual(len(matches), 0)

        # 2 reflections can match
        matches = dm.match_molecular_graphs(
            mol4, mol5, match_atoms_on=[], match_bonds_on=["bond_type"]
        )
        self.assertEqual(len(matches), 2)

        # 0 matches with hydrogens
        matches = dm.match_molecular_graphs(
            mol1, mol6, match_atoms_on=[], match_bonds_on=["bond_type"]
        )
        self.assertEqual(len(matches), 0)

        # 768 matches with hydrogens when specifying explicit_hs.
        # 768 = 6 rotations * 2 reflections * (2 reflection per carbon with Hs = 2^6)

        matches = dm.match_molecular_graphs(
            dm.add_hs(mol1), mol6, match_atoms_on=[], match_bonds_on=["bond_type"]
        )
        self.assertEqual(len(matches), 768)

    def test_mol(self):
        smiles1 = (
            "O=C(CC[C@@H]1NC(=O)N(CCc2c[nH]c3ccccc23)C1=O)N1C[C@@H]2C[C@H](C1)[C]1C=C[CH]C(=O)N1C2"
        )
        mol1 = dm.to_mol(smiles1, ordered=False)
        mol1_ordered = dm.to_mol(smiles1, ordered=True)

        smiles2 = "O=C1C2=C(N=CN2C)N(C(=O)N1C)C"
        mol2 = dm.to_mol(smiles2, ordered=False)
        mol2_ordered = dm.to_mol(smiles2, ordered=True)

        # Single match. Check if re-ordered atomic number match
        matches = dm.match_molecular_graphs(
            mol1, mol1_ordered, match_atoms_on=["atomic_num"], match_bonds_on=["bond_type"]
        )
        self.assertEqual(len(matches), 1)
        match = matches[0]
        atoms1 = [atom.GetAtomicNum() for atom in mol1.GetAtoms()]
        atoms1_ordered = [atom.GetAtomicNum() for atom in mol1_ordered.GetAtoms()]
        atoms1_re_ordered = [atoms1_ordered[match[ii]] for ii in range(len(match))]
        self.assertListEqual(atoms1, atoms1_re_ordered)

        # Single match. Check if re-ordered atomic number match
        matches = dm.match_molecular_graphs(
            mol2, mol2_ordered, match_atoms_on=["atomic_num"], match_bonds_on=["bond_type"]
        )
        self.assertEqual(len(matches), 1)
        match = matches[0]
        atoms2 = [atom.GetAtomicNum() for atom in mol2.GetAtoms()]
        atoms2_ordered = [atom.GetAtomicNum() for atom in mol2_ordered.GetAtoms()]
        atoms2_re_ordered = [atoms2_ordered[match[ii]] for ii in range(len(match))]
        self.assertListEqual(atoms2, atoms2_re_ordered)

        # Molecules don't match
        matches = dm.match_molecular_graphs(
            mol1, mol2, match_atoms_on=["atomic_num"], match_bonds_on=["bond_type"]
        )
        self.assertEqual(len(matches), 0)


class Test_reorder_mol_from_template(ut.TestCase):
    def test_reorder_mol_from_template(self):
        # No re-ordering because too many matches
        mol1 = dm.to_mol("C1CCCCC1", ordered=False)
        mol2 = dm.to_mol("C1=CC=CC=C1", ordered=True)
        mol_reordered = dm.reorder_mol_from_template(mol1, mol2)
        self.assertIsNone(mol_reordered)

        # Re-ordering with option ambiguous_match_mode
        mol1 = dm.to_mol("C1CCCCC1", ordered=False)
        mol2 = dm.to_mol("C1=CC=CC=C1", ordered=True)
        mol_reordered = dm.reorder_mol_from_template(mol1, mol2, ambiguous_match_mode="first")
        self.assertIsInstance(mol_reordered, dm.Mol)
        self.assertTrue(dm.same_mol(mol_reordered, mol1))

        # Re-ordering with option allow_ambiguous_hs_only should work since only hydrogens are ambiguous
        mol1 = dm.add_hs(dm.to_mol("C[N]1C=NC2=C1C(=O)N(C(=O)N2C)C", ordered=False))
        mol2 = dm.add_hs(dm.to_mol("C[N]1C=NC2=C1C(=O)N(C(=O)N2C)C", ordered=True))
        mol_reordered = dm.reorder_mol_from_template(mol1, mol2, ambiguous_match_mode="hs-only")
        self.assertIsInstance(mol_reordered, dm.Mol)
        self.assertTrue(dm.same_mol(mol_reordered, mol1))

        # Re-ordering with option allow_ambiguous_hs_only should not work since number of hydrogens are ambiguous
        mol1 = dm.add_hs(dm.to_mol("C1CCCCC1", ordered=False))
        mol2 = dm.add_hs(dm.to_mol("C1CCCCC1", ordered=True))
        mol_reordered = dm.reorder_mol_from_template(mol1, mol2, ambiguous_match_mode="hs-only")
        self.assertIsNone(mol_reordered)

        # Re-ordering with option allow_ambiguous_hs_only should not work since some non-hydrogens are ambiguous
        mol1 = dm.add_hs(dm.to_mol("C1CCCCC1", ordered=False))
        mol2 = dm.add_hs(dm.to_mol("C1=CC=CC=C1", ordered=True))
        mol_reordered = dm.reorder_mol_from_template(mol1, mol2, ambiguous_match_mode="hs-only")
        self.assertIsNone(mol_reordered)

        # Re-ordering with option ambiguous_match_mode=="best" should work
        mol1 = dm.to_mol("O=C(C)OC1=C(C=C(C=C1)CC)C(=O)O", ordered=False)
        mol2 = dm.to_mol("O=C(C)OC1=C(C=C(C=C1)C=C)C(=O)O", ordered=True)
        mol_reordered = dm.reorder_mol_from_template(
            mol1,
            mol2,
            enforce_atomic_num=True,
            enforce_bond_type=False,
            ambiguous_match_mode="best",
        )
        self.assertIsInstance(mol_reordered, dm.Mol)

        # Re-ordering with option ambiguous_match_mode=="best" should not work due to hydrogens
        mol1 = dm.to_mol("O=C(C)OC1=C(C=C(C=C1)CC)C(=O)O", ordered=False)
        mol2 = dm.to_mol("O=C(C)OC1=C(C=C(C=C1)C=C)C(=O)O", ordered=True)
        mol1, mol2 = dm.add_hs(mol1), dm.add_hs(mol2)
        mol_reordered = dm.reorder_mol_from_template(
            mol1,
            mol2,
            enforce_atomic_num=True,
            enforce_bond_type=False,
            ambiguous_match_mode="best",
        )
        self.assertIsNone(mol_reordered)

        # Re-ordering works with best, followed by first should work
        mol1 = dm.to_mol("O=C(C)OC1=C(C=C(C=C1)CC)C(=O)O", ordered=False)
        mol2 = dm.to_mol("O=C(C)OC1=C(C=C(C=C1)C=C)C(=O)O", ordered=True)
        mol_reordered = dm.reorder_mol_from_template(
            mol1,
            mol2,
            enforce_atomic_num=True,
            enforce_bond_type=False,
            ambiguous_match_mode="best-first",
        )
        self.assertIsInstance(mol_reordered, dm.Mol)
        self.assertEqual(
            sum(
                [
                    bond1.GetBondType() != bond2.GetBondType()
                    for bond1, bond2 in zip(mol_reordered.GetBonds(), mol2.GetBonds())
                ]
            ),
            1,
        )
        self.assertEqual(
            sum(
                [
                    atom1.GetAtomicNum() != atom2.GetAtomicNum()
                    for atom1, atom2 in zip(mol_reordered.GetAtoms(), mol2.GetAtoms())
                ]
            ),
            0,
        )

        # No-reordering because no matches
        mol1 = dm.to_mol("C1CCCCC1", ordered=False)
        mol2 = dm.to_mol("C1=CC=CC=C1", ordered=True)
        mol_reordered = dm.reorder_mol_from_template(mol1, mol2, enforce_bond_type=True)
        self.assertIsNone(mol_reordered)

        smiles3 = "O=C(CC[C@@H]3NC(=O)N(CCC1=C[N]([H])C2=C1C=CC=C2)C3=O)N4C[C@@H]6C[C@H](C4)C5C=CC([H])C(=O)N5C6"
        smiles3_variation = "O=C(CC[C@@H]3NC(=O)N(CCC1=C[N]([H])C2=C1C[N]CC2)C3=O)N4C[C@@H]6C[C@H](C4)C5C=CC([H])C(=O)N5C6"
        mol3 = dm.to_mol(smiles3, ordered=False)
        mol3_ordered = dm.to_mol(smiles3, ordered=True)
        mol3_variation = dm.to_mol(smiles3_variation, ordered=True)

        # Check re-ordering of molecules without enforcing
        mol3_reordered = dm.reorder_mol_from_template(
            mol3, mol3_ordered, enforce_atomic_num=True, enforce_bond_type=False
        )
        self.assertIsInstance(mol3_reordered, dm.Mol)
        atoms3 = [atom.GetAtomicNum() for atom in mol3.GetAtoms()]
        atoms3_ordered = [atom.GetAtomicNum() for atom in mol3_ordered.GetAtoms()]
        atoms3_reordered = [atom.GetAtomicNum() for atom in mol3_reordered.GetAtoms()]
        self.assertListEqual(atoms3_ordered, atoms3_reordered)
        self.assertFalse(all([atoms3[ii] == atoms3_ordered[ii] for ii in range(len(atoms3))]))

        # Check re-ordering of molecules without enforcing
        mol3_reordered = dm.reorder_mol_from_template(
            mol3, mol3_variation, enforce_atomic_num=False, enforce_bond_type=False
        )
        self.assertIsInstance(mol3_reordered, dm.Mol)
        atoms3 = [atom.GetAtomicNum() for atom in mol3.GetAtoms()]
        atoms3_ordered = [atom.GetAtomicNum() for atom in mol3_variation.GetAtoms()]
        atoms3_reordered = [atom.GetAtomicNum() for atom in mol3_reordered.GetAtoms()]
        equal_atoms3 = [atoms3_reordered[ii] == atoms3_ordered[ii] for ii in range(len(atoms3))]
        equal_atoms3_wrong = [atoms3_reordered[ii] == atoms3[ii] for ii in range(len(atoms3))]
        self.assertEqual(len(atoms3) - sum(equal_atoms3), 1)
        self.assertGreater(len(atoms3) - sum(equal_atoms3_wrong), 1)

        # Check re-ordering of molecules that fail due to enforcing
        mol3_reordered = dm.reorder_mol_from_template(
            mol3, mol3_variation, enforce_atomic_num=True, enforce_bond_type=False
        )
        self.assertIsNone(mol3_reordered)
        mol3_reordered = dm.reorder_mol_from_template(
            mol3, mol3_variation, enforce_atomic_num=False, enforce_bond_type=True
        )
        self.assertIsNone(mol3_reordered)


if __name__ == "__main__":
    ut.main()
