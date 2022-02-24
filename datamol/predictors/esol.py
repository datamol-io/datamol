from typing import Union

import pandas as pd


from .. import Mol

from ..descriptors.descriptors import clogp
from ..descriptors.descriptors import mw
from ..descriptors.descriptors import n_rotatable_bonds
from ..descriptors.descriptors import n_aromatic_atoms_proportion


_ESOL_INTERCEPT = 0.26121066137801696
_ESOL_COEF = {
    "mw": -0.0066138847738667125,
    "clogp": -0.7416739523408995,
    "n_rotatable_bonds": 0.003451545565957996,
    "n_aromatic_atoms_proportion": -0.42624840441316975,
}


def esol(mol: Mol):
    """Compute the solubility descriptor ESOL.

    Note that the intermediate descriptors will be computed on-the-fly. If you prefer
    precomputing those then you can use `esol_from_data`.

    Source: https://github.com/PatWalters/solubility/blob/d1536c58afe5e0e7ac4c96e2ffef496d5b98664b/esol.py
    """

    esol = (
        _ESOL_INTERCEPT
        + _ESOL_COEF["clogp"] * clogp(mol)
        + _ESOL_COEF["mw"] * mw(mol)
        + _ESOL_COEF["n_rotatable_bonds"] * n_rotatable_bonds(mol)
        + _ESOL_COEF["n_aromatic_atoms_proportion"] * n_aromatic_atoms_proportion(mol)
    )

    return esol


def esol_from_data(data: Union[pd.Series, pd.DataFrame, dict]):
    """Compute the solubility descriptor ESOL.

    `data` must contains the following intermediate descriptors:

    - `clogp`: `dm.descriptors.clogp`
    - `mw`: `dm.descriptors.mw`
    - `n_rotatable_bonds`: `dm.descriptors.n_rotatable_bonds`
    - `n_aromatic_atoms_proportion`: `dm.descriptors.n_aromatic_atoms_proportion`

    Source: https://github.com/PatWalters/solubility/blob/d1536c58afe5e0e7ac4c96e2ffef496d5b98664b/esol.py

    Args:
        data: A dataframe or series containing the intermediate descriptors.
    """

    esol = (
        _ESOL_INTERCEPT
        + _ESOL_COEF["clogp"] * data["clogp"]
        + _ESOL_COEF["mw"] * data["mw"]
        + _ESOL_COEF["n_rotatable_bonds"] * data["n_rotatable_bonds"]
        + _ESOL_COEF["n_aromatic_atoms_proportion"] * data["n_aromatic_atoms_proportion"]
    )

    return esol
