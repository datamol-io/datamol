**Added:**

* Add `dm.utils.fs.copy_dir()` to recursively copy directories across filesystems + tests.
* Add `dm.utils.fs.mkdir` + tests.
* Add a new `dm.descriptors` module with `compute_many_descriptors` and `batch_compute_many_descriptors` + tests.

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
