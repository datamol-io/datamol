**Added:**

* Add `dm.utils.fs.copy_dir()` to recursively copy directories across filesystems + tests.
* Add `dm.utils.fs.mkdir` + tests.
* Add a new `dm.descriptors` module with `compute_many_descriptors` and `batch_compute_many_descriptors` + tests.
* Add `dm.viz.match_substructure` to highlight one or more substructures in a list of molecules + tests. Note that the current function does not show different colors per match and submatch because of a limitation in `MolsToGridImage`. We plan to address this in a future version of datamol.
* Add a new `mcs` module backed by `rdkit.Chem.rdFMCS` with `find_mcs_with_details` and `find_mcs` functions + tests.

**Changed:**

* Set `fsspec` minimum version to `>=2021.09`.
* Pimp up `dm.utils.to_image` to make it more robust (don't fail on certain molecules due to incorrect aromaticity) and also propagate more drawing options to RDKit such as `legend_fontsize` and others.

**Deprecated:**

* <news item>

**Removed:**

* <news item>

**Fixed:**

* <news item>

**Security:**

* <news item>
