# `datamol`

`datamol` is a python library to work molecules. It's a layer built on top of [`rdkit`](https://www.rdkit.org/) and aims to be as light as possible.

- 🐍 Simple pythonic API
- ⚗️ Rdkit first: all you manipulate are `rdkit.Chem.Mol` objects.
- ✅ Manipulating molecules often rely on many options, datamol provides good defaults by design.
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

## Available modules

- `dm`: top-level module that contains common functions.
- `dm.actions`: functions to edit molecules.
- `dm.conformers`: generate and perform computation on conformers.
- `dm.data`: get some common data (mainly for dev purposes).
- `dm.fragment`: fragment molecules in a list of fragments.
- `dm.reactions`: functions to work with reactions.
- `dm.scaffold`: get representative scaffolds from a list of molecules.
- `dm.viz`: 2D/3D visualization functions.

## Install

Use conda:

```bash
conda install -c invivoai datamol
```

## Examples

See examples provided as a serie of [notebooks](./notebooks):

- [The Basics](docs/examples/The_Basics.ipynb)
- [Preprocessing Molecules](docs/examples/Preprocessing_Molecules.ipynb)
- [Cluster Molecules](docs/examples/Cluster_Molecules.ipynb)
- [Fragment and scaffold](docs/examples/Fragment_and_Scaffold.ipynb)

## Changelogs

See the latest changelogs at [CHANGELOG.rst](./CHANGELOG.rst).

## License

Under BSD license. See [LICENSE](LICENSE).

## Authors

See [AUTHORS.rst](./AUTHORS.rst).

## Development Lifecycle

## Build the documentation

You can build and serve the documentation locally with:

```bash
mkdocs serve
```

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
