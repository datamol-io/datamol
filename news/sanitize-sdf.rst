**Added:**

* Add a sanitize flag to `from_df`.
* Automatically detect the mol column in `from_df`.

**Changed:**

* Allow input a single molecule to `dm.to_sdf` instead of a list of mol.
* Preserve mol properties and the frist conformer in `dm.sanitize_mol`.
* Display a warning message when input mol has multiple conformers in `dm.sanitize_mol`.

**Deprecated:**

* <news item>

**Removed:**

* <news item>

**Fixed:**

* Remove call to `sanitize_mol` in `read_sdf`, instead use `sanitize=True` from RDKit.
* Remove the `mol` column from the mol properties in `from_df`. It also fixes `to_sdf`.

**Security:**

* <news item>
