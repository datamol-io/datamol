# Contribute

The below documents the development lifecycle of `datamol`.

## Setup a dev environment

```bash
conda create -n datamol
conda activate datamol

mamba env update -f env.yml

conda deactivate && conda activate datamol
pip install -e .
```

## Continuous Integration

`datamol` uses Github Actions to:

- **Build and test** `datamol`.
  - Multiple combinations of OS, Python and RDKit versions are tested.
- **Check** the code:
  - Formatting with `black`.
  - Static type check with `mypy`.
- **Documentation**: build and deploy the documentation on `master` and for every new git tag.

## Run tests

```bash
pytest
```

## Build the documentation

You can build and serve the documentation locally with:

```bash
# Regenerate the API doc
python -m datamol._mkdocs

# Build and serve the doc
mike serve
```

## Release a new version

- Install [rever](https://regro.github.io/rever-docs): `conda install -y rever`.
- Run check: `rever check`.
- Bump and release new version: `rever VERSION_NUMBER`.
- Releasing a new version will do the following things in that order:
  - Update `AUTHORS.rst`.
  - Update `CHANGELOG.rst`.
  - Bump the version number in `setup.py` and `_version.py`.
  - Add a git tag.
  - Push the git tag.
  - Add a new release on the GH repo associated with the git tag.
  - Update the appropriate feedstock to build a new conda package.
