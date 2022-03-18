**Added:**

* A new `dm.align` module with various functions to align a list of molecules. Use `dm.align.template_align` to align a molecule to a template and `dm.align.auto_align_many` to automatically partition and align a list of molecules.
* New descriptors: `formal_charge`
* New descriptors: `refractivity`
* New descriptors: `n_rigid_bonds`
* New descriptors: `n_stereo_centers`
* New descriptors: `n_charged_atoms`
* Add `dm.clear_props` to clear all the properties of a mol.
* Add a new dataset in addition to freesolv based on RDKit CDK2 at `dm.cdk2()`.
* Add `dm.strip_mol_to_core` to remove all R groups from a molecule.
* Add `dm.UNSPECIFIED_BOND`
* `dm.compute_ring_system` to extract the ring systems from a molecule.

**Changed:**

* Improve typing.
* Improve relative imports coverage.
* Adapt `dm.to_image` to use the `align` module.

**Deprecated:**

* <news item>

**Removed:**

* Remove a lot of `# type: ignore` as those can be error prone (hopefully the tests are here!)

**Fixed:**

* <news item>

**Security:**

* <news item>
