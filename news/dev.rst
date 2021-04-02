**Added:**

* `dm.copy_mol`
* `dm.set_mol_props`
* `dm.copy_mol_props`
* `dm.conformers.get_coords`
* `dm.conformers.center_of_mass`
* `dm.conformers.translate`
* `dm.enumerate_stereoisomers`
* `dm.enumerate_tautomers`
* `dm.atom_indices_to_mol`

**Changed:**

* rdkit fp to numpy array conversion is purely numpy-based now (x4 faster).
* Cleaning of various docstrings (removing explicit types).
* Clean various types.
* Allow `dm.to_image` instead of `dm.viz.to_image`
* Add atom indices drawing option to `dm.to_image`
* Allow to smiles to fail (default is to not fail but return None as before).
* Add CXSmiles bool flag to to_smiles.
* Rename utils.paths to utils.fs

**Deprecated:**

* <news item>

**Removed:**

* <news item>

**Fixed:**

* Scaffold tests for new rdkit version
* Conformer cluster tests for new rdkit version

**Security:**

* <news item>
