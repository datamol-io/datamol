from typing import Callable
from typing import Dict
from typing import List
from typing import Union
from typing import Optional

import functools

import pandas as pd

from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

from .. import Mol
from ..utils.jobs import parallelized

from .descriptors import mw
from .descriptors import fsp3
from .descriptors import n_lipinski_hba
from .descriptors import n_lipinski_hbd
from .descriptors import n_rings
from .descriptors import n_hetero_atoms
from .descriptors import n_heavy_atoms
from .descriptors import n_rotatable_bonds
from .descriptors import n_radical_electrons
from .descriptors import tpsa
from .descriptors import qed
from .descriptors import clogp
from .descriptors import sas
from .descriptors import n_aliphatic_carbocycles
from .descriptors import n_aliphatic_heterocyles
from .descriptors import n_aliphatic_rings
from .descriptors import n_aromatic_carbocycles
from .descriptors import n_aromatic_heterocyles
from .descriptors import n_aromatic_rings
from .descriptors import n_saturated_carbocycles
from .descriptors import n_saturated_heterocyles
from .descriptors import n_saturated_rings


def any_rdkit_descriptor(name: str) -> Callable:
    """Return a descriptor function by name either from
    `rdkit.Chem import Descriptors` or `rdkit.Chem.rdMolDescriptors`.

    Args:
        name: Descriptor name.
    """
    fn = getattr(Descriptors, name, None)

    if fn is None:
        fn = getattr(rdMolDescriptors, name, None)

    if fn is None:
        raise ValueError(f"Descriptor {name} not found.")

    return fn


_DEFAULT_PROPERTIES_FN = {
    "mw": mw,
    "fsp3": fsp3,
    "n_lipinski_hba": n_lipinski_hba,
    "n_lipinski_hbd": n_lipinski_hbd,
    "n_rings": n_rings,
    "n_hetero_atoms": n_hetero_atoms,
    "n_heavy_atoms": n_heavy_atoms,
    "n_rotatable_bonds": n_rotatable_bonds,
    "n_radical_electrons": n_radical_electrons,
    "tpsa": tpsa,
    "qed": qed,
    "clogp": clogp,
    "sas": sas,
    "n_aliphatic_carbocycles": n_aliphatic_carbocycles,
    "n_aliphatic_heterocyles": n_aliphatic_heterocyles,
    "n_aliphatic_rings": n_aliphatic_rings,
    "n_aromatic_carbocycles": n_aromatic_carbocycles,
    "n_aromatic_heterocyles": n_aromatic_heterocyles,
    "n_aromatic_rings": n_aromatic_rings,
    "n_saturated_carbocycles": n_saturated_carbocycles,
    "n_saturated_heterocyles": n_saturated_heterocyles,
    "n_saturated_rings": n_saturated_rings,
}


def compute_many_descriptors(
    mol: Mol,
    properties_fn: Optional[Dict[str, Union[Callable, str]]] = None,
    add_properties: bool = True,
) -> dict:
    """Compute a list of opiniated molecular properties.

    Args:
        mol: A molecule.
        properties_fn: A list of functions that compute properties. If None,
            a default list of properties is used. If the function is a string,
            `dm.descriptors.any_descriptor()` is used to retrieve the descriptor
            function.
        add_properties: Whether to add the computed properties to the default list.

    Returns:
        Computed properties as a dict.
    """

    if properties_fn is None:
        properties_fn = _DEFAULT_PROPERTIES_FN
    elif add_properties:
        [properties_fn.setdefault(k, v) for k, v in _DEFAULT_PROPERTIES_FN.items()]

    props = {}
    for k, v in properties_fn.items():

        if isinstance(v, str):
            v = any_rdkit_descriptor(v)

        props[k] = v(mol)

    return props


def batch_compute_many_descriptors(
    mols: List[Mol],
    properties_fn: Optional[Dict[str, Union[Callable, str]]] = None,
    add_properties: bool = True,
    n_jobs: int = 1,
    batch_size: Optional[int] = None,
    progress: bool = False,
    progress_leave: bool = True,
) -> pd.DataFrame:
    """Compute a list of opiniated molecular properties on a list of molecules.

    Args:
        mols: A list of molecules.
        properties_fn: A list of functions that compute properties. If None,
            a default list of properties is used. If the function is a string,
            `dm.descriptors.any_descriptor()` is used to retrieve the descriptor
            function.
        add_properties: Whether to add the computed properties to the default list.

    Returns:
        A dataframe of computed properties with one row per input molecules.
    """

    compute_fn = functools.partial(
        compute_many_descriptors,
        properties_fn=properties_fn,
        add_properties=add_properties,
    )

    props = parallelized(
        compute_fn,
        mols,
        batch_size=batch_size,
        progress=progress,
        n_jobs=n_jobs,
        tqdm_kwargs=dict(leave=progress_leave),
    )
    return pd.DataFrame(props)
