**Added:**

* Add `dm.utils.fs.copy_dir()` to recursively copy directories across filesystems + tests.
* Add `dm.utils.fs.mkdir` + tests.
* Add a new `dm.descriptors` module with `compute_many_descriptors` and `batch_compute_many_descriptors` + tests.
* Add `dm.viz.match_substructure` to highlight one or more substructures in a list of molecules + tests. Note that the current function does not show different colors per match and submatch because of a limitation in `MolsToGridImage`. We plan to address this in a future version of datamol.
* Add a new `mcs` module backed by `rdkit.Chem.rdFMCS` with `find_mcs_with_details` and `find_mcs` functions + tests.
* Add a new function `dm.viz.utils.align_2d_coordinates` to align 2d coordinates of molecules using either a given pattern or MCS.
* Add `dm.canonical_tautomer` to canonicalize tautomers.
* Add `dm.remove_stereochemistry()`.

**Changed:**

* Set `fsspec` minimum version to `>=2021.9`.
* Pimp up `dm.utils.to_image` to make it more robust (don't fail on certain molecules due to incorrect aromaticity) and also propagate more drawing options to RDKit such as `legend_fontsize` and others.
* Add a new `align` argument in `dm.to_image()` to align the 2d coordinates of the molecules.

**Deprecated:**

* <news item>

**Removed:**

* <news item>

**Fixed:**

* <news item>

**Security:**

* <news item>
