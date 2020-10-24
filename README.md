# `datamol`

**IMPORTANT: Still WIP. DO NOT USE.**

`datamol` is a pythonic library to manipulate molecules. It's a layer built on top of [`rdkit`](https://www.rdkit.org/) and aims to be as light as possible.

- Simple pythonic API
- Rdkit first: all you manipulate are `rdkit.Chem.Mol` objects.
- Performance matters: built-in efficient parallelization when possible with optional progress bar.
- Modern IO: out-of-the-box support for remote paths using `fsspec` to read and write multiple formats (sdf, xlsx, csv, etc).

## datamol's API

_NOTE(hadim): the below snippet is important. It allows people to have a quick idea of datamol's API. Let's make it nice (and short)!._

```python
import datamol as dm

# Convenient functions
mol = dm.to_mol("O=C(C)Oc1ccccc1C(=O)O", sanitize=True)
fp = dm.as_fp(mol)

mol_with_conformers = dm.add_confs(mol)

# Easy IO
mols = dm.read_sdf("s3://my-awesome-data-lake/smiles.sdf")

df = dm.read_csv("/home/data/dataset.csv")
```

## Changelogs

See the latest changelogs at [CHANGELOG.rst](./CHANGELOG.rst).

## Install

Use conda:

```bash
conda install -c invivoai datamol
```

## Documentation

TODO (try [mkdocs?](https://www.mkdocs.org/))

## Examples

See examples provided as a serie of [notebooks](./notebooks) (TODO).

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
