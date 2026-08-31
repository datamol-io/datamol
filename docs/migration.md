# Migrating to Datamol 1.x

Datamol 1.x refreshes the supported scientific Python stack without adding new
features. The goal is to preserve Datamol's public behaviour while removing
compatibility code for dependency versions that are no longer maintained.

## Supported runtime

- Python 3.11 or newer is required.
- RDKit 2024.09 or newer is required.
- NumPy 1.26, pandas 2.2 and scikit-learn 1.4 are the minimum supported series.

Install the core package with `python -m pip install datamol`. Optional dependency
groups are available for I/O (`datamol[io]`) and visualisation (`datamol[viz]`).

## Dependency-driven output changes

Recent RDKit releases may choose a different, but chemically equivalent, canonical
SMILES or CXSMILES representation. Conformer embedding and clustering can also
retain a different number of conformers between RDKit releases. Code should compare
chemical structure and conformer invariants rather than version-specific strings,
energies or cluster counts.

pandas 3 removed the `verbose` argument from `pandas.read_csv`. For compatibility,
`datamol.open_df(..., verbose=...)` accepts and ignores that argument with a
deprecation warning when the installed pandas version no longer supports it. Other
unknown reader arguments continue to raise an error.

## Development and releases

Create the development environment with `mamba env create -n datamol -f env.yml`.
The CI tests the supported Python and RDKit series, tutorial notebooks, documentation,
formatting and package distributions separately. Releases are built from published
GitHub releases and uploaded to PyPI with short-lived OpenID Connect credentials.
