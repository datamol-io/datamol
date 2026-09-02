# Contribute

The below documents the development lifecycle of Datamol.

## Setup a dev environment

```bash
uv sync --all-extras
```

This creates an isolated `.venv` and installs Datamol in editable mode with the
test, documentation, I/O, visualisation and development extras. `env.yml`
remains available for contributors who need a Conda environment.

Run the same formatting and lint checks used by CI before committing:

```bash
uv run pre-commit run --all-files
```

## Setup a dev environment with dev container

This repository is setup to use [dev container](https://docs.github.com/en/codespaces/setting-up-your-project-for-codespaces/introduction-to-dev-containers). You can use it locally with VSCode or any editor supporting dev containers as well as on GitHub Codespaces.

The env is based on the Micromamba Docker image.

## Continuous Integration

Datamol uses Github Actions to:

- **Build and test** `datamol`.
  - Python 3.11 through 3.14 and the supported RDKit release series are tested.
  - The current stack is also tested on Linux x86-64, Windows x86-64,
    macOS Apple Silicon and macOS Intel.
  - Tutorial notebooks run in a dedicated job so failures are easier to diagnose.
- **Check** the code:
  - Formatting with `black`.
  - Linting with `ruff`.
  - Building and validating the wheel and source distribution.
- **Documentation**: build on pull requests, and deploy from `main` and successful manual releases.

## Run tests

```bash
uv run python -m pytest -m "not integration"
uv run python -m pytest -m integration --no-cov -n 0
```

The first command is the fast core suite. The second executes the maintained
tutorials and is the same integration command used by the dedicated GitHub
Actions job.

## Build the documentation

You can build and serve the documentation locally with:

```bash
# Build and serve the doc
uv run mike serve
```

### Multi-versioning

The documentation is built for every pull request and deployed for pushes to
`main` and successful manual releases using [mike](https://github.com/jimporter/mike).
Everything is automated using GitHub Actions.

## Release a new version

Run the `release` action manually from `main`, with the intended version and
`dry-run` unchecked when ready to publish. The existing `PYPI_API_TOKEN`
secret authenticates PyPI uploads. Tests, package validation and documentation
must pass first. See the [release guide](releasing.md) for the rehearsal,
prerelease and recovery steps.

The existing conda-forge feedstock remains the Conda release channel. After a
PyPI release, conda-forge's update bot proposes the new version; maintainers
review its dependencies and tests in the feedstock rather than duplicating a
Conda upload inside this repository's release workflow.
