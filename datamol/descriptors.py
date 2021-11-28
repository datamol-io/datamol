from typing import Callable
from typing import Dict
from typing import List

import sys
import functools
import os

import pandas as pd

from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import RDConfig

sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
import sascorer  # type:ignore

from . import Mol
from .utils.jobs import parallelized

mw = rdMolDescriptors.CalcExactMolWt
fsp3 = rdMolDescriptors.CalcFractionCSP3
n_hba = rdMolDescriptors.CalcNumHBA
n_hbd = rdMolDescriptors.CalcNumHBD
n_rings = rdMolDescriptors.CalcNumRings
n_hetero_atoms = rdMolDescriptors.CalcNumHeteroatoms
n_heavy_atoms = Descriptors.HeavyAtomCount  # type: ignore
n_rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds
n_aliphatic_rings = Descriptors.NumAliphaticRings  # type: ignore
n_aromatic_rings = Descriptors.NumAromaticRings  # type: ignore
n_saturated_rings = Descriptors.NumSaturatedRings  # type: ignore
n_H_acceptors = Descriptors.NumHAcceptors  # type: ignore
n_H_donors = Descriptors.NumHDonors  # type: ignore
n_radical_electrons = Descriptors.NumRadicalElectrons
tpsa = rdMolDescriptors.CalcTPSA
qed = Descriptors.qed
clogp = Descriptors.MolLogP  # type: ignore
sas = sascorer.calculateScore

_DEFAULT_PROPERTIES_FN = {
    "mw": mw,
    "fsp3": fsp3,
    "n_hba": n_hba,
    "n_hbd": n_hbd,
    "n_rings": n_rings,
    "n_hetero_atoms": n_hetero_atoms,
    "n_heavy_atoms": n_heavy_atoms,
    "n_rotatable_bonds": n_rotatable_bonds,
    "n_aliphatic_rings": n_aliphatic_rings,
    "n_aromatic_rings": n_aromatic_rings,
    "n_saturated_rings": n_saturated_rings,
    "n_H_acceptors": n_H_acceptors,
    "n_H_donors": n_H_donors,
    "n_radical_electrons": n_radical_electrons,
    "tpsa": tpsa,
    "qed": qed,
    "clogp": clogp,
    "sas": sas,
}


def compute_many_descriptors(
    mol: Mol,
    properties_fn: Dict[str, Callable] = None,
    add_properties: bool = True,
) -> dict:
    """Compute a list of opiniated molecular properties.

    Args:
        mol: A molecule.
        properties_fn: A list of functions that compute properties. If None,
            a default list of properties is used.
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
        props[k] = v(mol)

    return props


def batch_compute_many_descriptors(
    mols: List[Mol],
    properties_fn: Dict[str, Callable] = None,
    add_properties: bool = True,
    n_jobs: int = 1,
    batch_size: int = None,
    progress: bool = False,
    progress_leave: bool = True,
) -> pd.DataFrame:
    """Compute a list of opiniated molecular properties on a list of molecules.

    Args:
        mols: A list of molecules.
        properties_fn: A list of functions that compute properties. If None,
            a default list of properties is used.
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
