==================
Datamol Changelogs
==================

.. current developments

v0.4.5
====================

**Added:**

* Add missing appdirs dependency
* Add missing appdirs dependency

**Fixed:**

* Propagate tqdm_kwargs for parallel (was only done for sequential)

**Authors:**

* Hadrien Mary



v0.4.4
====================

**Added:**

* Add `tqdm_kwargs` to `dm.utils.JobRunner()`
* Add `tqdm_kwargs` to `dm.utils.parallelized()`

**Changed:**

* Propagate `job_kwargs` to dm.utils.parallelized()`

**Authors:**

* Hadrien Mary



v0.4.3
====================

**Added:**

* Add a DOI so datamol can get properly cited.
* Better doc about compat and CI
* Add a datamol Mol type: `dm.Mol` identical to `Chem.rdchem.Mol`

**Changed:**

* Bump test coverage from 70% to 80%.

**Authors:**

* DeepSource Bot
* Hadrien Mary
* deepsource-autofix[bot]



v0.4.2
====================

**Added:**

* More tests for the `dm.similarity` modules + check against RDKit equivalent methods.
* `dm.same_mol(mol1, mol2)` to check whether 2 molecules are the same based on their InChiKey.

**Changed:**

* use `scipy` in `dm.similarity.pdist()`.
* Raise an error when a molecule is invalid in `dm.similarity.pdist/cdist`.

**Deprecated:**

* `dm.similarity.pdist()` nows returns only the dist matrix without the `valid_idx` vector.

**Fixed:**

* A bug returning an inconsistent dist matrix with `dm.similarity.pdist()`.

**Authors:**

* Hadrien Mary



v0.4.1
====================

**Changed:**

* A better and manually curated API documentation.

**Authors:**

* Hadrien Mary



v0.4.0
====================

**Added:**

* Add support for more fingerprint types.
* Two utility functions for molar concentration conversion: `dm.molar_to_log()` and `dm.log_to_molar()`.
* Add the `dm.utils.fs` module to work with any type of paths (remote or local).

**Authors:**

* Hadrien Mary



v0.3.9
====================

**Added:**

* Add a sanitize flag to `from_df`.
* Automatically detect the mol column in `from_df`.
* Add `add_hs` arg to `sanitize_mol`.

**Changed:**

* Allow input a single molecule to `dm.to_sdf` instead of a list of mol.
* Preserve mol properties and the frist conformer in `dm.sanitize_mol`.
* Display a warning message when input mol has multiple conformers in `dm.sanitize_mol`.

**Fixed:**

* Remove call to `sanitize_mol` in `read_sdf`, instead use `sanitize=True` from RDKit.
* Remove the `mol` column from the mol properties in `from_df`. It also fixes `to_sdf`.

**Authors:**

* Hadrien Mary



v0.3.8
====================

**Changed:**

* Propagate `sanitize` and `strict_parsing` to `dm.read_sdf`.

**Authors:**

* Hadrien Mary
* Ishan Kumar
* michelml



v0.3.7
====================

**Fixed:**

* Fix again and hopefully the last time google analytics.

**Authors:**

* Hadrien Mary



v0.3.6
====================

**Changed:**

* Add s3fs and gcsfs as hard dep

**Authors:**

* Hadrien Mary



v0.3.5
====================

**Authors:**

* Hadrien Mary
* michelml



v0.3.4
====================

**Authors:**

* Hadrien Mary



v0.3.3
====================

**Changed:**

* New logo.

**Authors:**

* Hadrien Mary



v0.3.2
====================

**Fixed:**

* Fixed typo in readme

**Authors:**

* Emmanuel Noutahi
* Hadrien Mary



v0.3.1
====================

**Authors:**

* Hadrien Mary



v0.3.0
====================

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
* Integrate pandatools into `dm.to_df`.
* Build a mol column from smiles in read_csv and read_excel
* Rename `dm.sanitize_best` to `dm.sanitize_first`
*

**Fixed:**

* Scaffold tests for new rdkit version
* Conformer cluster tests for new rdkit version

**Authors:**

* Hadrien Mary
* Therence1
* michelml
* mike



v0.2.12
====================

**Fixed:**

* Tqdm progress bar update on completion of job and not submission

**Authors:**

* Emmanuel Noutahi



v0.2.11
====================

**Changed:**

* Make ipywidgets an optional dep.

**Authors:**

* Hadrien Mary



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




