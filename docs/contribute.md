# Contribute

The below documents the development lifecycle of Datamol.

## Setup a dev environment

```bash
mamba env create -n datamol -f env.yml
mamba activate datamol
```

The environment installs Datamol in editable mode with the test, documentation,
I/O, visualisation and development dependency groups.

Run the same formatting and lint checks used by CI before committing:

```bash
pre-commit run --all-files
```

## Setup a dev environment with dev container

This repository is setup to use [dev container](https://docs.github.com/en/codespaces/setting-up-your-project-for-codespaces/introduction-to-dev-containers). You can use it locally with VSCode or any editor supporting dev containers as well as on GitHub Codespaces.

The env is based on the Micromamba Docker image.

## Continuous Integration

Datamol uses Github Actions to:

- **Build and test** `datamol`.
  - Python 3.11 through 3.14 and the supported RDKit release series are tested.
  - The current stack is also tested on Linux, macOS and Windows.
  - Tutorial notebooks run in a dedicated job so failures are easier to diagnose.
- **Check** the code:
  - Formatting with `black`.
  - Linting with `ruff`.
  - Building and validating the wheel and source distribution.
- **Documentation**: build on pull requests, and deploy from `main` and published releases.

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

### Multi-versioning

The documentation is built for every pull request and deployed for pushes to
`main` and published releases using [mike](https://github.com/jimporter/mike).
Everything is automated using GitHub Actions.

## Release a new version

Create and publish a GitHub release from a version tag after the release commit has
landed on `main`. The [`release` workflow](https://github.com/datamol-io/datamol/actions/workflows/release.yml)
builds and validates both distributions, tests the wheel in a clean environment,
and publishes it through PyPI Trusted Publishing. The `pypi` GitHub environment
and matching PyPI trusted publisher must be configured once by a repository
administrator; no long-lived PyPI API token is stored in GitHub.
