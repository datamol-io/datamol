**Added:**

* Support for `max_num_mols` in `dm.read_sdf()`. Useful when files are large and debugging code.
* Support for returning the invalid molecules in `dm.read_sdf`. Useful when we need to know which one failed.
* Support for more compression formats when reading SDF files using `fssep.open(..., compression="infer")`.
* Add `CODEOWNERS` file.
* Add `dm.descriptors.n_spiro_atoms` and `dm.descriptors.n_stereo_centers_unspecified`.

**Changed:**

* Overload output types for `dm.read_sdf` and `dm.data.*`.

**Deprecated:**

* <news item>

**Removed:**

* <news item>

**Fixed:**

* <news item>

**Security:**

* <news item>
