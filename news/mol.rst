**Added:**

* Add multiple utilities to work with mapped SMILES with hydrogens.
* Add `dm.clear_atom_props()` to remove atom's properties.
* Add `dm.clear_atom_map_number()` to remove the atom map number property.
* Add `dm.get_atom_positions()` to retrieve the atomic positions of a conformer of a molecule.
* Add `dm.set_atom_positions()` to add a new confomer to a molecule given a list of atomic positions.

**Changed:**

* Add new arguments to `dm.to_mol`: `allow_cxsmiles`, `parse_name`, `remove_hs` and `strict_cxsmiles`. Refers to the docstring for the details.
* Set `copy` to `True` by default to `dm.atom_indices_to_mol()`.
* Allow to specify the property keys to clear in `dm.clear_mol_props()`. If not set, the original default beahviour is to clear everything.

**Deprecated:**

* <news item>

**Removed:**

* <news item>

**Fixed:**

* <news item>

**Security:**

* <news item>
