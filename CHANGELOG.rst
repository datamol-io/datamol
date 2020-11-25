=====================
pantagruel Change Log
=====================

.. current developments

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




