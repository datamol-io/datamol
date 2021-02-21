<img src="docs/images/logo-title-200.png" height="80px">

# The rdkit-based molecular Python library

<div align="center">

[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/datamol)](https://pypi.org/project/datamol/)
[![PyPI Status](https://badge.fury.io/py/datamol.svg)](https://badge.fury.io/py/datamol)
[![Conda](https://img.shields.io/conda/v/conda-forge/datamol?label=conda&color=success)](https://anaconda.org/conda-forge/datamol)
[![ReadTheDocs](https://readthedocs.org/projects/datamol/badge/?version=stable)](https://datamol.readthedocs.io/en/stable/)
[![license](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://github.com/datamol-org/datamol/blob/master/LICENSE)

</div>

`datamol` is a python library to work with molecules. It's a layer built on top of [RDKit](https://www.rdkit.org/) and aims to be as light as possible.

- 🐍 Simple pythonic API
- ⚗️ Rdkit first: all you manipulate are `rdkit.Chem.Mol` objects.
- ✅ Manipulating molecules often rely on many options, datamol provides good defaults by design.
- 🧠 Performance matters: built-in efficient parallelization when possible with optional progress bar.
- 🕹️ Modern IO: out-of-the-box support for remote paths using `fsspec` to read and write multiple formats (sdf, xlsx, csv, etc).

## Entrypoints

- Website: https://datamol.io
- Documentation: https://datamol.readthedocs.io/en/stable/

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
mol = dm.standardized_mol(mol)

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

## Install

Use conda:

```bash
mamba install -c invivoai datamol
```

## Changelogs

See the latest changelogs at [CHANGELOG.rst](./CHANGELOG.rst).

## License

Under the Apache-2.0 license. See [LICENSE](LICENSE).

## Authors

See [AUTHORS.rst](./AUTHORS.rst).
