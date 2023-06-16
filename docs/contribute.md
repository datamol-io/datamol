# Contribute

The below documents the development lifecycle of Datamol.

## Setup a dev environment

```bash
mamba env create -n datamol -f env.yml
mamba activate datamol
pip install -e .
```

## Setup a dev environment with dev container

This repository is setup to use [dev container](https://docs.github.com/en/codespaces/setting-up-your-project-for-codespaces/introduction-to-dev-containers). You can use it locally with VSCode or any editor supporting dev containers as well as on GitHub Codespaces.

The env is based on the Micromamba Docker image.

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

The process is fully automated by executing the [`release` GH Action](https://github.com/datamol-io/datamol/actions/workflows/release.yml).
