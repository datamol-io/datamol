# Migrating to Datamol 1.x

Datamol 1.x refreshes the supported scientific Python stack without adding new
features. The goal is to preserve Datamol's public behaviour while removing
compatibility code for dependency versions that are no longer maintained.

## Supported runtime

- Python 3.11 or newer is required.
- RDKit 2024.09 or newer is required.
- NumPy 1.26, pandas 2.2 and scikit-learn 1.4 are the minimum supported series.

Install the core package with `python -m pip install datamol`. Optional dependency
groups are available for I/O (`datamol[io]`), visualisation (`datamol[viz]`),
and SELFIES conversion (`datamol[selfies]`). Matplotlib, Pillow, nglview and
SELFIES no longer make every core installation heavier. HTTP I/O remains
available in the default installation.

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

`fold_count_fp` now accepts both sparse and explicit RDKit bit vectors, as its
public type contract already indicated. Datamol remains the canonical home for
general fingerprint folding and conformer alignment; downstream packages such
as Molfeat reuse these primitives rather than maintaining divergent copies.

`template_align(..., auto_select_coord_gen=True)` now scopes RDKit's CoordGen
preference to the alignment call. It no longer changes the process-wide RDKit
depiction preference for subsequent calls.

## Development and releases

Create the development environment with `uv sync --all-extras`; `env.yml`
remains a supported Conda alternative. The CI uses the same uv-based install
path and tests the supported Python and RDKit series on Linux x86-64, Windows
x86-64, macOS Apple Silicon and macOS Intel. Tutorial notebooks,
documentation, formatting and package distributions run separately. Releases
are built from published GitHub releases, smoke-tested as both wheels and source
distributions, attested, and uploaded to PyPI with short-lived OpenID Connect
credentials. The conda-forge feedstock continues to consume PyPI releases via
its update bot and remains a supported downstream distribution.
