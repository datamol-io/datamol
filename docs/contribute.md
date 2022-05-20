# Contribute

The below documents the development lifecycle of Datamol.

## Setup a dev environment

```bash
mamba create -n datamol -f env.yml
conda activate datamol
pip install -e .
```

## Continuous Integration

Datamol uses Github Actions to:

- **Build and test** `datamol`.
  - Multiple combinations of OS, Python and RDKit versions are tested.
- **Check** the code:
  - Formatting with `black`.
  - Static type check with `mypy`.
- **Documentation**: build and deploy the documentation on `main` and for every new git tag.

## Run tests

```bash
pytest
```

## Build the documentation

You can build and serve the documentation locally with:

```bash
# Build and serve the doc
mike serve
```

### Multi-versionning

The doc is built for eash push on `main` and every git tags using [mike](https://github.com/jimporter/mike). Everything is automated using Github Actions. Please refer to the official mike's documentation for the details.

## Release a new version

- Run check: `rever check`.
- Bump and release new version: `rever VERSION_NUMBER`.
- Releasing a new version will do the following things in that order:
  - Update `AUTHORS.rst`.
  - Update `CHANGELOG.rst`.
  - Bump the version number in `setup.py` and `_version.py`.
  - Add a git tag.
  - Push the git tag.
  - Add a new release on the GH repo associated with the git tag.
