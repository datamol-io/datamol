# Datamol Changelogs

## v0.10.3

**Authors:**

- Hadrien Mary

## v0.10.2

**Authors:**

- Hadrien Mary

## v0.10.1

**Changed:**

- Different docs style, with Tabs on the top for `Overview`, `Usage`,
  `Tutorials`, `API`, `Contribute`, `License`

**Authors:**

- DomInvivo

## v0.10.0

**Added:**

- Add two new fucntions `dm.open_df` and `dm.save_df` that
  automatically open and save a dataframe to a file. This is a
  convenience function that automatically handles the file format
  based on the file extension. The functions are also able to handle
  compressed files.

**Changed:**

- All the datamol modules and objects are now lazy loaded. It means
  that loading now happens on-demand. Preliminary tests suggest the
  datamol import time decreases by 20-fold (from 1s to 50ms on a
  regular Ubuntu laptop) without affecting the subsequent calls to the
  modules and objects. This is a major improvement for the datamol
  usability. This new behaviour is enabled by default but can be
  disabled by setting the environment variable
  `DATAMOL_DISABLE_LAZY_LOADING` to `1`.
- Move the fs module to its dedicated section in the docs. Fix #160.

**Removed:**

- Remove unused, broken and uncovered
  `datamol.fragment.assemble_fragment_iter()` function.

**Authors:**

- Hadrien Mary
- dessygil

## v0.9.8

**Authors:**

- Hadrien Mary

## v0.9.7

**Authors:**

- Hadrien Mary

## v0.9.4

**Fixed:**

- Fix wrong image output for lasso viz function. Make it consistent
  with `dm.to_image()` and rdkit.
- Avoid global `IPython` import so it's not an hard datamol
  dependency.
- Add `importlib-resources` dep in the datamol pypi package.

**Authors:**

- Hadrien Mary

## v0.9.3

**Added:**

- added a feature that highlights substructures of 2D molecular images

**Changed:**

- Update CNAME to docs.datamol.io
- Replace all occurrences of doc.datamol.io by docs.datamol.io
- Switch from `pkg_resources` to `importlib.resources` for loading
  resources.
- Enable python 3.11 on the CI.
- Relocatem `datamol/data.py` to `datamol/data/\_\_init\_\_.py`.

**Fixed:**

- Color bug of the search input bar

**Authors:**

- Emmanuel Noutahi
- Hadrien Mary
- Honoré Hounwanou
- dessygil

## v0.9.2

**Added:**

- A multi-mol2 file reader that converts into rdkit objects

**Fixed:**

- Updated the logging in `\_sanifix4.py` to use the RDKit logger

**Authors:**

- Cas
- Hadrien Mary
- Pakman450
- Therence1

## v0.9.1

**Changed:**

- moved `CODE_OF_CONDUCT.md`, `CODEOWNDERS`, `CONTRIBUTING.md` and
  `SECURITY.md` to `.github/` dir
- Improve and automate the release process.
- Adapt the logo and colors to the new branding.
- Replace `datamol-org` to `datamol-io` everywhere in the codebase due
  to GH org rename.

**Authors:**

- Hadrien Mary
- Saurav Maheshkar

## v0.9.0

**Changed:**

- Add `TypeAlias` types to `datamol.types.\*`.
- Drop `setup.py` in favour of `pyproject.toml` only.
- Replace unmaintained `appdirs` by maintained `platformdirs`.
- Enable weekly tests on `main` branch.

**Fixed:**

- Add missing fcfp func in fingerprint functions dict

**Authors:**

- Hadrien Mary
- michelml

## v0.8.8

**Added:**

- Add PDB read/writer functions: `dm.to_pdbblock()`,
  `dm.read_pdbblock()`, `dm.read_pdbfile()`, `dm.to_pdbfile()`

**Changed:**

- Improve output type in `to_df`.\`

**Authors:**

- Hadrien Mary

## v0.8.7

**Added:**

- Add multiple utilities to work with mapped SMILES with hydrogens.
- Add `dm.clear_atom_props()` to remove atom's properties.
- Add `dm.clear_atom_map_number()` to remove the atom map number
  property.
- Add `dm.get_atom_positions()` to retrieve the atomic positions of a
  conformer of a molecule.
- Add `dm.set_atom_positions()` to add a new confomer to a molecule
  given a list of atomic positions.

**Changed:**

- Add new arguments to \`dm.to_mol\`: `allow_cxsmiles`, `parse_name`,
  `remove_hs` and `strict_cxsmiles`. Refers to the docstring for the
  details.
- Set `copy` to `True` by default to `dm.atom_indices_to_mol()`.
- Allow to specify the property keys to clear in
  `dm.clear_mol_props()`. If not set, the original default beahviour
  is to clear everything.

**Authors:**

- Hadrien Mary

## v0.8.6

**Fixed:**

- Ensure rdkit 2021.03 works with latest datamol. The support is not
  "official" but only a single function must be adapted so it's ok.

**Authors:**

- Hadrien Mary

## v0.8.5

**Added:**

- Support for `max_num_mols` in `dm.read_sdf()`. Useful when files are
  large and debugging code.
- Support for returning the invalid molecules in `dm.read_sdf`. Useful
  when we need to know which one failed.
- Support for more compression formats when reading SDF files using
  `fssep.open(..., compression="infer")`.
- Add `CODEOWNERS` file.
- Add `dm.descriptors.n_spiro_atoms` and
  `dm.descriptors.n_stereo_centers_unspecified`.

**Changed:**

- Overload output types for `dm.read_sdf` and `dm.data.\*`.
- Reduce tests duration (especially in CI).

**Authors:**

- DomInvivo
- Hadrien Mary

## v0.8.4

**Changed:**

- Add a comment recommending to not use the SMI file format.

**Fixed:**

- Fix a bug when reading a remote file with `dm.read_smi()`.

**Authors:**

- Hadrien Mary

## v0.8.3

**Added:**

- Parallelization to `to_df` for faster conversion to dataframe

**Fixed:**

- Error in docs

**Authors:**

- Emmanuel Noutahi

## v0.8.2

**Fixed:**

- Fix a typo in a tutorial.

**Authors:**

- Hadrien Mary
- Valence-JonnyHsu

## v0.8.1

**Changed:**

- Remove the `rdkit` dependency in the setup.py to prevent pip to
  always override the conda rdkit package. See
  <https://github.com/rdkit/rdkit/pull/2690#issuecomment-1295375416>
  for context.

**Authors:**

- Hadrien Mary

## v0.8.0

**Added:**

- `dm.Atom` and `dm.Bond` types.
- Add RDKit as a pypi dep.
- Add `datamol.hash_mol()` based on `rdkit.Chem.RegistrationHash`.

**Changed:**

- RDKit 2022.09: use `Draw.shouldKekulize` instead of
  `Draw.\_okToKekulizeMol`.
- RDKit 2022.09: don't use `dm.convert.\_ChangeMoleculeRendering` for
  RDKit \>=2022.09.

**Authors:**

- Hadrien Mary

## v0.7.18

**Added:**

- Added argument product_index in `select_reaction_output`. It allows
  to return all products and a product of interest by the index.
- Updated unit tests.

**Authors:**

- Lu Zhu

## v0.7.17

**Added:**

- Added a new chemical reaction module for rdkit chemical reactions
  and attachment manipulations.

**Fixed:**

**Authors:**

- Hadrien Mary
- Lu Zhu

## v0.7.16

**Changed:**

- Bump upstream GH actions versions.
- `dm.fs.copy_dir` now uses the internal fsspec `copy` when the two
  source and destination fs are the same. It makes the copy much
  faster.

**Fixed:**

- Use `os.PathLike` to recognize a broader range of string-based path
  inputs in the `dm.fs` module. It prevents file objects such as
  `py.\_path.local.LocalPath` not being recognized as path.

**Authors:**

- Hadrien Mary

## v0.7.15

**Fixed:**

- Missing header in the fragment tutorial.

**Authors:**

- Hadrien Mary
- Valence-JonnyHsu

## v0.7.14

**Added:**

- Add `with_atom_indices` to `dm.to_smiles`. If enable, atom indices
  will be added to the SMILES.

**Changed:**

- Changed the default for `dm.fs.is_file()` from
  ` True`` to \ `False\`.
- Refactor the API doc to breakdown all the submodules in individual
  doc. Thanks to @MichelML for the suggestion.
- Re-enable pipy activity in rever.

**Fixed:**

- Minor typo in the documentation of `dm.conformers.generate()`

**Authors:**

- Cas
- Hadrien Mary
- Valence-JonnyHsu

## v0.7.13

**Added:**

- New aligning tutorials.

**Removed:**

- `rdkit` dep from pypi (the dep is only on the conda forge package)

**Fixed:**

- Grammar in tutorials.

**Authors:**

- Hadrien Mary
- Valence-JonnyHsu

## v0.7.12

**Fixed:**

- Fix minor typos in tutorials

**Authors:**

- Hadrien Mary
- michelml

## v0.7.11

**Added:**

- Add configurations for dev containers based on the micromamba Docker
  image. More informations about dev container at
  <https://docs.github.com/en/codespaces/setting-up-your-project-for-codespaces/introduction-to-dev-containers>.
- support for two additional forcefields: MMFF94s with and without
  electrostatic component
- energies output along with delta-energy to lowest energy conformer

**Changed:**

- API of dm.conformers.generate() to support choice of forcefield. In
  addition ewindow and eratio flags added to reject high energy
  conformers, either on absoute scale, or as ratio to rotatable bonds
- Revamped all the datamol tutorials and add new tutorials. Huge
  thanks to @Valence-jonnyhsu for leading the refactoring of the
  datamol tutorials.
- Improve documentation for `dm.standardize_mol()`
- Multiple various docstring and typing improvments.
- Embed the cdk2.sdf and [solubility]()\*.sdf files within the datamol
  package to prevent issue with the RDKit config dir.
- Enable strict mode on the documentation to prevent any issues and
  inconsistency with the types and docstrings of datamol.
- Refactor micromamba CI to use latest and simplify it.

**Removed:**

- Remove unused and unmaintained `dm.actions` and `dm.reactions`
  module.
- Remove `copy` args from `add_hs` and `remove_hs` (RDKit already
  returns copies).

**Fixed:**

- Errors in ECFP fingerprints that computes FCFP instead of ECFP.

**Authors:**

- Emmanuel Noutahi
- Hadrien Mary
- Matt

## v0.7.10

**Added:**

- New possibilities for ambiguous matching of molecules in the
  function `reorder_mol_from_template`

**Changed:**

- Replaced `allow_ambiguous_hs_only` by the option `"hs_only"` for the
  `ambiguous_match_mode` parameter
- `ambiguous_match_mode` is now a String, no longer a bool.

**Deprecated:**

- `allow_ambiguous_hs_only` is no longer deprecated, but without
  warning since the feature is brand new.
- Same for `ambiguous_match_mode` being a bool.

**Authors:**

- DomInvivo
- Hadrien Mary

## v0.7.9

**Added:**

- `datamol.graph.match_molecular_graphs`, with unit-tests
- `datamol.graph.reorder_mol_from_template`, with unit-tests

**Changed:**

- Typing in `datamol.graph.py`, changed `rdkit.Chem.rdchem.Mol` to
  `dm.Mol`

**Deprecated:**

- NOTHING

**Removed:**

- NOTHING

**Fixed:**

- NOTHING

**Security:**

- NOTHING

**Authors:**

- DomInvivo
- Emmanuel Noutahi

## v0.7.8

**Fixed:**

- Bug in `dm.conformer.generate()` when multiple conformers had equal
  energies
- Fix the documentation.

**Authors:**

- Cas
- Hadrien Mary

## v0.7.7

**Added:**

- Add `dm.read_molblock()` and `dm.to_molblock()` functions.
- Add `dm.to_xlsx()` function.

**Fixed:**

- Fix the API doc.

**Authors:**

- Hadrien Mary

## v0.7.6

**Changed:**

- Add `joblib_batch_size` in `dm.parallelized_with_batches()` to be
  able to control the joblib batch size (which is different than the
  `dm.parallelized_with_batches` batch size.
- Various small improvements for unit tests.

**Authors:**

- Hadrien Mary

## v0.7.5

**Added:**

- Add `dm.parallelized_with_batches()` to parallelize workload with a
  function that take a batch of inputs.

**Authors:**

- Hadrien Mary

## v0.7.4

**Changed:**

- Don't import `sasscorer` by default but only during the call to
  `dm.descriptors.sas(mol)`

**Authors:**

- Hadrien Mary

## v0.7.3

**Changed:**

- Use micromamba during CI.
- Add CI tests for RDKit=2022.03.
- Adapt a test to new rdkit version.

**Fixed:**

- typing for what is returned by dm.align.template_align

**Authors:**

- Hadrien Mary
- michelml

## v0.7.2

**Changed:**

- allow_r_groups option in dm.align.auto_align_many

**Removed:**

- should_align

**Authors:**

- Hadrien Mary
- michelml

## v0.7.1

**Added:**

- A new `dm.align` module with various functions to align a list of
  molecules. Use `dm.align.template_align` to align a molecule to a
  template and `dm.align.auto_align_many` to automatically partition
  and align a list of molecules.
- New descriptors: `formal_charge`
- New descriptors: `refractivity`
- New descriptors: `n_rigid_bonds`
- New descriptors: `n_stereo_centers`
- New descriptors: `n_charged_atoms`
- Add `dm.clear_props` to clear all the properties of a mol.
- Add a new dataset in addition to freesolv based on RDKit CDK2 at
  `dm.cdk2()`.
- Add `dm.strip_mol_to_core` to remove all R groups from a molecule.
- Add `dm.UNSPECIFIED_BOND`
- `dm.compute_ring_system` to extract the ring systems from a
  molecule.

**Changed:**

- Improve typing.
- Improve relative imports coverage.
- Adapt `dm.to_image` to use the `align` module.

**Removed:**

- Remove a lot of `\# type: ignore` as those can be error prone
  (hopefully the tests are here!)

**Authors:**

- Hadrien Mary

## v0.7.0

**Added:**

- Add `dm.conformers.keep_conformers` in order to only keep one or
  multiple conformers from a molecules.

**Changed:**

- Change the conformer generation arguments to use
  `useRandomCoords=True` by default.
- Start using explicit `Optional` instead of implicit `Optional` for
  typing.
- Start using relative imports instead of absolute ones.
- When conformers are not minimized, sort them by energy (can be
  turned to False).

**Removed:**

- Remove `fallback_to_random_coords` argument from
  `generate_conformers`.

**Authors:**

- Hadrien Mary

## v0.6.9

**Added:**

- Support for selfies\<2.0.0 in tests

**Changed:**

- Behaviour of all _inchi_ functions to return None with a warning
  instead of silently returning an empty string
- Order of str evaluation on convertion function. `isinstance(str)` is
  now evaluated before `is None`

**Fixed:**

- Bug in unique_id making this evaluation falling back on
  'd41d8cd98f00b204e9800998ecf8427e' on unsupported inputs. Instead
  None is returned now

**Authors:**

- Emmanuel Noutahi

## v0.6.8

**Changed:**

- Add `remove_hs` flag in `dm.read_sdf()`.

**Authors:**

- Hadrien Mary

## v0.6.7

**Added:**

- Add `dm.descriptors.n_aromatic_atoms`
- Add `dm.descriptors.n_aromatic_atoms_proportion`
- Add `dm.predictors.esol`
- Add `dm.predictors.esol_from_data`

**Changed:**

- Make `descriptors` a folder (backward compatible).
- Rename `any_descriptor` to `any_rdkit_descriptor` to be more
  explicit.

**Authors:**

- Hadrien Mary

## v0.6.6

**Added:**

- Add `dm.conformers.align_conformers()` to align the conformers of a
  list of molecules.

**Changed:**

- New lower bound rdkit version to `\>=2021.09`. See #81 for details.

**Authors:**

- Hadrien Mary

## v0.6.5

**Fixed:**

- Catch too long integer values in `set_mol_props` and switch to
  `SetDoubleProp` instead of `SetIntProp`

**Authors:**

- Hadrien Mary

## v0.6.4

**Changed:**

- Expose the clean_it flag when enumerating stereoisomers.

**Authors:**

- Hadrien Mary
- Julien Horwood

## v0.6.3

**Added:**

- Parameters allowing to customize or ignore failures when running the
  conformer generation.

**Changed:**

- When the conformer embedding fails, it will now optionally fall back
  to using random coordinates.

**Authors:**

- Hadrien Mary
- Julien Horwood

## v0.6.2

**Added:**

- Add a new `total` arg in `dm.parallelized()` (only useful when the
  `progress` is set to `True`)

**Changed:**

- Prevent ` tqdm_kwargs`` collision in \ `dm.parallelized()\`.

**Authors:**

- Hadrien Mary

## v0.6.1

**Added:**

- Add `dm.to_inchi_non_standard()` and `dm.to_inchikey_non_standard()`
  in order to generate InChi values that are sensitive to tautomerism
  as well as undefined stereoisomery.
- Add `dm.unique_id` to generate unique molecule identifiers based on
  `dm.to_inchikey_non_standard`

**Changed:**

- Add `use_non_standard_inchikey` flag argument to `dm.same_mol`.

**Authors:**

- Hadrien Mary

## v0.6.0

**Added:**

- Add `dm.utils.fs.copy_dir()` to recursively copy directories across
  filesystems + tests.
- Add `dm.utils.fs.mkdir` + tests.
- Add a new `dm.descriptors` module with `compute_many_descriptors`
  and `batch_compute_many_descriptors` + tests.
- Add `dm.viz.match_substructure` to highlight one or more
  substructures in a list of molecules + tests. Note that the current
  function does not show different colors per match and submatch
  because of a limitation in `MolsToGridImage`. We plan to address
  this in a future version of datamol.
- Add a new `mcs` module backed by `rdkit.Chem.rdFMCS` with `find_mcs`
  function + tests.
- Add a new function `dm.viz.utils.align_2d_coordinates` to align 2d
  coordinates of molecules using either a given pattern or MCS.
- Add `dm.canonical_tautomer` to canonicalize tautomers.
- Add `dm.remove_stereochemistry()`.
- Add a `bond_line_width` arg to `to_image`.
- Add `dm.atom_list_to_bond()`
- Add `enable` flag to `dm.without_rdkit_log()`
- Add a tutorial about the filesystem module.
- Add a tutorial about the viz module (still incomplete).
- Add `dm.substructure_matching_bonds` to perform a standard
  substructure match but also return the matching bonds instead of
  only the matching atoms.
- Add new `dm.isomers` module + move relevant functions from `dm.mol`
  to `dm.isomers`
- Add `dm.add_hs` and `dm.remove` to add and remove hydrogens from
  molecules.

**Changed:**

- Set `fsspec` minimum version to `\>=2021.9`.
- Pimp up `dm.utils.to_image` to make it more robust (don't fail on
  certain molecules due to incorrect aromaticity) and also propagate
  more drawing options to RDKit such as `legend_fontsize` and others.
- Add a new `align` argument in `dm.to_image()` to align the 2d
  coordinates of the molecules.
- In `dm.to_image`, `use_svg` is now set to `True` by default.
- Change the default `mol_size` from 200 to 300 in `to_image`.
- Link `datamol.utils.fs` to `datamol.fs`.
- Change default `chunk_size` in `copy_file` from 2048 to 1024 \* 1024
  (1MB).
- Support parallel chunked distances computation in
  `dm.similarity.cdist`

**Authors:**

- Hadrien Mary

## v0.5.0

**Changed:**

- The default git branch is now `main`
- `appdirs` is now an hard dep.
- Change CI to use rdkit `\[2021.03, 2021.09\]` and add the info the
  readme and doc.

**Fixed:**

- Test related to SELFIES to make it work with the latest 2.0 version.
- `dm.to_mol` accept `mol` as input but the specified type was only
  `str`.

**Authors:**

- Hadrien Mary

## v0.4.11

**Fixed:**

- Force the input value(s) of `dm.molar.log_to_molar` to be a float
  since power of integers are not allowed.

**Authors:**

- Hadrien Mary

## v0.4.10

**Removed:**

- `py.typed` file that seems unused beside confusing static analyzer
  tools.

**Authors:**

- Hadrien Mary

## v0.4.9

**Added:**

- `to_smarts` for exporting molecule objects as SMARTS
- `from_smarts` for reading molecule from SMARTS string

**Changed:**

- Allow exporting smiles in kekule representaiton
- `to_smarts` is properly renamed into `smiles_as_smarts`

**Authors:**

- Emmanuel Noutahi

## v0.4.8

**Removed:**

- Revert batch_size fix to use default joblib instead

**Fixed:**

- Issue #58: sequence bug in parallel.

**Authors:**

- Emmanuel Noutahi

## v0.4.7

**Added:**

- Add a new function to measure execution time
  `dm.utils.perf.watch_duration`.

**Changed:**

- Add a `batch_size` option to `dm.utils.parallelized`. The default
  behaviour `batch_size=None` is unchanged and so 100% backward
  compatible.

**Authors:**

- Hadrien Mary

## v0.4.6

**Changed:**

- `get_protocol` is more general

**Fixed:**

- Bug in fs.glob due to protocol being a list

**Authors:**

- Emmanuel Noutahi

## v0.4.5

**Added:**

- Add missing appdirs dependency
- Add missing appdirs dependency

**Fixed:**

- Propagate tqdm_kwargs for parallel (was only done for sequential)

**Authors:**

- Hadrien Mary

## v0.4.4

**Added:**

- Add `tqdm_kwargs` to `dm.utils.JobRunner()`
- Add `tqdm_kwargs` to `dm.utils.parallelized()`

**Changed:**

- Propagate `job_kwargs` to dm.utils.parallelized()\`

**Authors:**

- Hadrien Mary

## v0.4.3

**Added:**

- Add a DOI so datamol can get properly cited.
- Better doc about compat and CI
- Add a datamol Mol type: `dm.Mol` identical to `Chem.rdchem.Mol`

**Changed:**

- Bump test coverage from 70% to 80%.

**Authors:**

- DeepSource Bot
- Hadrien Mary
- deepsource-autofix\[bot\]

## v0.4.2

**Added:**

- More tests for the `dm.similarity` modules + check against RDKit
  equivalent methods.
- `dm.same_mol(mol1, mol2)` to check whether 2 molecules are the same
  based on their InChiKey.

**Changed:**

- use `scipy` in `dm.similarity.pdist()`.
- Raise an error when a molecule is invalid in
  `dm.similarity.pdist/cdist`.

**Deprecated:**

- `dm.similarity.pdist()` nows returns only the dist matrix without
  the `valid_idx` vector.

**Fixed:**

- A bug returning an inconsistent dist matrix with
  `dm.similarity.pdist()`.

**Authors:**

- Hadrien Mary

## v0.4.1

**Changed:**

- A better and manually curated API documentation.

**Authors:**

- Hadrien Mary

## v0.4.0

**Added:**

- Add support for more fingerprint types.
- Two utility functions for molar concentration conversion:
  `dm.molar_to_log()` and `dm.log_to_molar()`.
- Add the `dm.utils.fs` module to work with any type of paths (remote
  or local).

**Authors:**

- Hadrien Mary

## v0.3.9

**Added:**

- Add a sanitize flag to `from_df`.
- Automatically detect the mol column in `from_df`.
- Add `add_hs` arg to `sanitize_mol`.

**Changed:**

- Allow input a single molecule to `dm.to_sdf` instead of a list of
  mol.
- Preserve mol properties and the frist conformer in
  `dm.sanitize_mol`.
- Display a warning message when input mol has multiple conformers in
  `dm.sanitize_mol`.

**Fixed:**

- Remove call to `sanitize_mol` in `read_sdf`, instead use
  `sanitize=True` from RDKit.
- Remove the `mol` column from the mol properties in `from_df`. It
  also fixes `to_sdf`.

**Authors:**

- Hadrien Mary

## v0.3.8

**Changed:**

- Propagate `sanitize` and `strict_parsing` to `dm.read_sdf`.

**Authors:**

- Hadrien Mary
- Ishan Kumar
- michelml

## v0.3.7

**Fixed:**

- Fix again and hopefully the last time google analytics.

**Authors:**

- Hadrien Mary

## v0.3.6

**Changed:**

- Add s3fs and gcsfs as hard dep

**Authors:**

- Hadrien Mary

## v0.3.5

**Authors:**

- Hadrien Mary
- michelml

## v0.3.4

**Authors:**

- Hadrien Mary

## v0.3.3

**Changed:**

- New logo.

**Authors:**

- Hadrien Mary

## v0.3.2

**Fixed:**

- Fixed typo in readme

**Authors:**

- Emmanuel Noutahi
- Hadrien Mary

## v0.3.1

**Authors:**

- Hadrien Mary

## v0.3.0

**Added:**

- `dm.copy_mol`
- `dm.set_mol_props`
- `dm.copy_mol_props`
- `dm.conformers.get_coords`
- `dm.conformers.center_of_mass`
- `dm.conformers.translate`
- `dm.enumerate_stereoisomers`
- `dm.enumerate_tautomers`
- `dm.atom_indices_to_mol`

**Changed:**

- rdkit fp to numpy array conversion is purely numpy-based now (x4
  faster).

- Cleaning of various docstrings (removing explicit types).

- Clean various types.

- Allow `dm.to_image` instead of `dm.viz.to_image`

- Add atom indices drawing option to `dm.to_image`

- Allow to smiles to fail (default is to not fail but return None as
  before).

- Add CXSmiles bool flag to to_smiles.

- Rename utils.paths to utils.fs

- Integrate pandatools into `dm.to_df`.

- Build a mol column from smiles in read_csv and read_excel

- Rename `dm.sanitize_best` to `dm.sanitize_first`

- **Fixed:**

- Scaffold tests for new rdkit version

- Conformer cluster tests for new rdkit version

**Authors:**

- Hadrien Mary
- Therence1
- michelml
- mike

## v0.2.12

**Fixed:**

- Tqdm progress bar update on completion of job and not submission

**Authors:**

- Emmanuel Noutahi

## v0.2.11

**Changed:**

- Make ipywidgets an optional dep.

**Authors:**

- Hadrien Mary

## v0.2.10

**Changed:**

- Propagate more options to dm.reorder_atoms.

**Authors:**

- Hadrien Mary

## v0.2.9

**Added:**

- `dm.pick_centroids` for picking a set of centroid molecules using
  various algorithm
- `dm.assign_to_centroids` for clustering molecules based on
  precomputed centroids.

**Changed:**

- Make `add_hs` optional in `conformers.generate` and removed them
  when `add_hs` is True. Explicit hydrogens will be lost.

**Fixed:**

- Doc string of `dm.pick_diverse`

**Authors:**

- Emmanuel Noutahi
- Hadrien Mary

## v0.2.8

**Added:**

- Added outfile to viz.to_image

**Changed:**

- Replace ete3 by networkx due to GPL licensing.
- Fix some typos in docs.

**Fixed:**

- Null pointer exception during conformers generation.

**Authors:**

- Emmanuel Noutahi
- Hadrien Mary
- Honoré Hounwanou
- michelml

## v0.2.7

**Added:**

- Add a test to monitor datamol import duration.

**Changed:**

- Add rms cutoff option during conformers generation.
- Refactor conformer cluster function.

**Authors:**

- Hadrien Mary

## v0.2.6

**Added:**

- Include stub files for rdkit generated using stubgen from mypy.

**Authors:**

- Hadrien Mary

## v0.2.5

**Added:**

- Add `to_smi` and `from_smi` in the IO module.
- Support filelike object in io module.
- Add kekulization to to_mol

**Changed:**

- Switch tests of the IO module to regular functions.

**Deprecated:**

- In the IO module, use `urlpath` instead of `file_uri` to follow
  `fsspec` conventions.

**Fixed:**

- Fix bug in read_excel where sheet_name wasnt being used.

**Authors:**

- Emmanuel Noutahi
- Hadrien Mary

## v0.2.4

**Changed:**

- Constraint rdkit to 2020.09 to get `rdBase.LogStatus()`

**Authors:**

- Hadrien Mary

## v0.2.3

**Changed:**

- Better rdkit log disable/enable.

**Authors:**

- Hadrien Mary

## v0.2.2

**Added:**

- Test that execute the notebooks.

**Fixed:**

- Force rdkit \>=2020.03.6 to avoid thread-related bug in
  `rdMolStandardize`

**Authors:**

- Hadrien Mary

## v0.2.1

**Added:**

- Add `cdist` function to compute tanimoto sim between two list of
  molecules.

**Fixed:**

- Fix a bug in `dm.from_df` when the dataframe has a size of zero.

**Authors:**

- Hadrien Mary

## v0.2.0

**Added:**

- Add all the common sanitize functions.
- Add the 2_Preprocessing_Molecules notebook.
- Add fragment module.
- Add scaffold module.
- Add cluster module.
- Add assemble module.
- Add actions module.
- Add reactions module.
- Add dm.viz.circle_grid function
- Add doc with mkdocs

**Authors:**

- Hadrien Mary

## v0.1.2

**Authors:**

- Hadrien Mary

## v0.1.1

**Authors:**

## v0.1.0

**Added:**

- first release!
