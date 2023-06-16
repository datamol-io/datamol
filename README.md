<div align="center">
    <img src="docs/images/logo-title.svg" width="100%">
</div>

<p align="center">
    <b>datamol - molecular processing made easy</b> <br />
</p>
<p align="center">
  <a href="https://docs.datamol.io/stable/" target="_blank">
      Docs
  </a> |
  <a href="https://datamol.io/" target="_blank">
      Homepage
  </a>
</p>

---

[![DOI](https://zenodo.org/badge/341603042.svg)](https://zenodo.org/badge/latestdoi/341603042)
[![Binder](http://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/datamol-io/datamol/main?urlpath=lab/tree/docs/tutorials/The_Basics.ipynb)
[![PyPI](https://img.shields.io/pypi/v/datamol)](https://pypi.org/project/datamol/)
[![Conda](https://img.shields.io/conda/v/conda-forge/datamol?label=conda&color=success)](https://anaconda.org/conda-forge/datamol)
[![PyPI - Downloads](https://img.shields.io/pypi/dm/datamol)](https://pypi.org/project/datamol/)
[![Conda](https://img.shields.io/conda/dn/conda-forge/datamol)](https://anaconda.org/conda-forge/datamol)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/datamol)](https://pypi.org/project/datamol/)
[![license](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://github.com/datamol-io/datamol/blob/main/LICENSE)
[![GitHub Repo stars](https://img.shields.io/github/stars/datamol-io/datamol)](https://github.com/datamol-io/datamol/stargazers)
[![GitHub Repo stars](https://img.shields.io/github/forks/datamol-io/datamol)](https://github.com/datamol-io/datamol/network/members)
[![Codecov](https://codecov.io/gh/datamol-io/datamol/branch/main/graph/badge.svg?token=2ETG8SA7IG)](https://codecov.io/gh/datamol-io/datamol)

Datamol is a python library to work with molecules. It's a layer built on top of [RDKit](https://www.rdkit.org/) and aims to be as light as possible.

- 🐍 Simple pythonic API
- ⚗️ RDKit first: all you manipulate are `rdkit.Chem.Mol` objects.
- ✅ Manipulating molecules often relies on many options; Datamol provides good defaults by design.
- 🧠 Performance matters: built-in efficient parallelization when possible with an optional progress bar.
- 🕹️ Modern IO: out-of-the-box support for remote paths using `fsspec` to read and write multiple formats (sdf, xlsx, csv, etc).

## Try Online

Visit [![Binder](http://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/datamol-io/datamol/main?urlpath=lab/tree/docs/tutorials/The_Basics.ipynb) and try Datamol online.

## Documentation

Visit <https://docs.datamol.io>.

## Installation

Use conda:

```bash
mamba install -c conda-forge datamol
```

## Quick API Tour

```python
import datamol as dm

# Common functions
mol = dm.to_mol("O=C(C)Oc1ccccc1C(=O)O", sanitize=True)
fp = dm.to_fp(mol)
selfies = dm.to_selfies(mol)
inchi = dm.to_inchi(mol)

# Standardize and sanitize
mol = dm.to_mol("O=C(C)Oc1ccccc1C(=O)O")
mol = dm.fix_mol(mol)
mol = dm.sanitize_mol(mol)
mol = dm.standardize_mol(mol)

# Dataframe manipulation
df = dm.data.freesolv()
mols = dm.from_df(df)

# 2D viz
legends = [dm.to_smiles(mol) for mol in mols[:10]]
dm.viz.to_image(mols[:10], legends=legends)

# Generate conformers
smiles = "O=C(C)Oc1ccccc1C(=O)O"
mol = dm.to_mol(smiles)
mol_with_conformers = dm.conformers.generate(mol)

# 3D viz (using nglview)
dm.viz.conformers(mol, n_confs=10)

# Compute SASA from conformers
sasa = dm.conformers.sasa(mol_with_conformers)

# Easy IO
mols = dm.read_sdf("s3://my-awesome-data-lake/smiles.sdf", as_df=False)
dm.to_sdf(mols, "gs://data-bucket/smiles.sdf")
```

## How to cite

Please cite Datamol if you use it in your research: [![DOI](https://zenodo.org/badge/341603042.svg)](https://zenodo.org/badge/latestdoi/341603042).

## Compatibilities

Version compatibilities are an essential topic for production-software stacks. We are cautious about documenting compatibility between `datamol`, `python` and `rdkit`.

See below the associated versions of Python and RDKit, for which a minor version of Datamol **has been tested** during its whole lifecycle. _It does not mean other combinations does not work but that those are not tested._

| `datamol` | `python`            | `rdkit`                       |
| --------- | ------------------- | ----------------------------- |
| `0.11.x`  | `[3.9, 3.10, 3.11]` | `[2022.09, 2023.03]`          |
| `0.10.x`  | `[3.9, 3.10, 3.11]` | `[2022.03, 2022.09]`          |
| `0.9.x`   | `[3.9, 3.10, 3.11]` | `[2022.03, 2022.09]`          |
| `0.8.x`   | `[3.8, 3.9, 3.10]`  | `[2021.09, 2022.03, 2022.09]` |
| `0.7.x`   | `[3.8, 3.9]`        | `[2021.09, 2022.03]`          |
| `0.6.x`   | `[3.8, 3.9]`        | `[2021.09]`                   |
| `0.5.x`   | `[3.8, 3.9]`        | `[2021.03, 2021.09]`          |
| `0.4.x`   | `[3.8, 3.9]`        | `[2020.09, 2021.03]`          |
| `0.3.x`   | `[3.8, 3.9]`        | `[2020.09, 2021.03]`          |

## CI Status

The CI runs tests and performs code quality checks for the following combinations:

- The three major platforms: Windows, OSX and Linux.
- The two latest Python versions.
- The two latest RDKit versions.

|                                         | `main`                                                                                                                                                                    |
| --------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Lib build & Testing                     | [![test](https://github.com/datamol-io/datamol/actions/workflows/test.yml/badge.svg)](https://github.com/datamol-io/datamol/actions/workflows/test.yml)                   |
| Code Sanity (linting and type analysis) | [![code-check](https://github.com/datamol-io/datamol/actions/workflows/code-check.yml/badge.svg)](https://github.com/datamol-io/datamol/actions/workflows/code-check.yml) |
| Documentation Build                     | [![doc](https://github.com/datamol-io/datamol/actions/workflows/doc.yml/badge.svg)](https://github.com/datamol-io/datamol/actions/workflows/doc.yml)                      |

## License

Under the Apache-2.0 license. See [LICENSE](LICENSE).
