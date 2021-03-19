=====================
pantagruel Change Log
=====================

.. current developments

v0.2.10
====================

**Changed:**

* Propagate more options to dm.reorder_atoms.

**Authors:**

* Hadrien Mary



v0.2.9
====================

**Added:**

* `dm.pick_centroids` for picking a set of centroid molecules using various algorithm
* `dm.assign_to_centroids` for clustering molecules based on precomputed centroids.

**Changed:**

* Make `add_hs` optional in `conformers.generate` and removed them when `add_hs` is True. Explicit hydrogens will be lost.

**Fixed:**

* Doc string of `dm.pick_diverse`

**Authors:**

* Emmanuel Noutahi
* Hadrien Mary



v0.2.8
====================

**Added:**

* Added outfile to viz.to_image

**Changed:**

* Replace ete3 by networkx due to GPL licensing.
* Fix some typos in docs.

**Fixed:**

* Null pointer exception during conformers generation.

**Authors:**

* Emmanuel Noutahi
* Hadrien Mary
* Honoré Hounwanou
* michelml



v0.2.7
====================

**Added:**

* Add a test to monitor datamol import duration.

**Changed:**

* Add rms cutoff option during conformers generation.
* Refactor conformer cluster function.

**Authors:**

* Hadrien Mary



v0.2.6
====================

**Added:**

* Include stub files for rdkit generated using stubgen from mypy.

**Authors:**

* Hadrien Mary



v0.2.5
====================

**Added:**

* Add `to_smi` and `from_smi` in the IO module.
* Support filelike object in io module.
* Add kekulization to to_mol

**Changed:**

* Switch tests of the IO module to regular functions.

**Deprecated:**

* In the IO module, use `urlpath` instead of `file_uri` to follow `fsspec` conventions.

**Fixed:**

* Fix bug in read_excel where sheet_name wasnt being used.

**Authors:**

* Emmanuel Noutahi
* Hadrien Mary



v0.2.4
====================

**Changed:**

* Constraint rdkit to 2020.09 to get `rdBase.LogStatus()`

**Authors:**

* Hadrien Mary



v0.2.3
====================

**Changed:**

* Better rdkit log disable/enable.

**Authors:**

* Hadrien Mary



v0.2.2
====================

**Added:**

* Test that execute the notebooks.

**Fixed:**

* Force rdkit >=2020.03.6 to avoid thread-related bug in `rdMolStandardize`

**Authors:**

* Hadrien Mary



v0.2.1
====================

**Added:**

* Add `cdist` function to compute tanimoto sim between two list of molecules.

**Fixed:**

* Fix a bug in `dm.from_df` when the dataframe has a size of zero.

**Authors:**

* Hadrien Mary



v0.2.0
====================

**Added:**

* Add all the common sanitize functions.
* Add the 2_Preprocessing_Molecules notebook.
* Add fragment module.
* Add scaffold module.
* Add cluster module.
* Add assemble module.
* Add actions module.
* Add reactions module.
* Add dm.viz.circle_grid function
* Add doc with mkdocs

**Authors:**

* Hadrien Mary



v0.1.2
====================

**Authors:**

* Hadrien Mary



v0.1.1
====================

**Authors:**




v0.1.0
====================

**Added:**

* first release!

**Authors:**




