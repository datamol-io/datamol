# `datamol`

`datamol` is a python library to work molecules. It's a layer built on top of [`rdkit`](https://www.rdkit.org/) and aims to be as light as possible.

- 🐍 Simple pythonic API
- ⚗️ Rdkit first: all you manipulate are `rdkit.Chem.Mol` objects.
- ✅ Good default: manipulating molecules often rely on a lot different options, datamol provides good default options.
- 🧠 Performance matters: built-in efficient parallelization when possible with optional progress bar.
- 🕹️ Modern IO: out-of-the-box support for remote paths using `fsspec` to read and write multiple formats (sdf, xlsx, csv, etc).

## API

```python
import datamol as dm

# Common functions
mol = dm.to_mol("O=C(C)Oc1ccccc1C(=O)O", sanitize=True)
fp = dm.to_fp(mol)
selfies = dm.to_selfies(mol)
inchi = dm.to_inchi(mol)

# Generate conformers
mol_with_conformers = dm.add_confs(mol)

# Easy IO
mols = dm.read_sdf("s3://my-awesome-data-lake/smiles.sdf", as_df=False)
dm.to_sdf(mols, "gs://data-bucket/smiles.sdf")
```

## Install

Use conda:

```bash
conda install -c invivoai datamol
```

## Examples

See examples provided as a serie of [notebooks](./notebooks):

1. [The Basics](notebooks/1_The_Basics.ipynb)

## Documentation

TODO (try [mkdocs?](https://www.mkdocs.org/))

## Changelogs

See the latest changelogs at [CHANGELOG.rst](./CHANGELOG.rst).

## License

Under BSD license. See [LICENSE](LICENSE).

## Authors

See [AUTHORS.rst](./AUTHORS.rst).

## Release a new version

- Install [rever](https://regro.github.io/rever-docs): `conda install -y rever`.
- Run check: `rever check`.
- Bump and release new version: `rever VERSION_NUMBER`.
- Releasing a new version will do the following things in that order:
  - Update [AUTHORS.rst](./AUTHORS.rst).
  - Update [CHANGELOG.rst](./CHANGELOG.rst).
  - Bump the version number in `setup.py` and `_version.py`.
  - Add a git tag.
  - Push the git tag.
  - Add a new release on the GH repo associated with the git tag.
  - Update the appropriate feedstock to build a new conda package.
