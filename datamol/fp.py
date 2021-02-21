from typing import Union
from typing import Optional
from typing import List

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import DataStructs

import numpy as np

import datamol as dm


def fp_to_array(fp: DataStructs.ExplicitBitVect, dtype: type = int) -> np.ndarray:
    """Convert rdkit fingerprint to numpy array."""
    if isinstance(fp, np.ndarray):
        return fp
    arr = np.zeros((0,), dtype=dtype)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr


def to_fp(
    mol: Union[str, Chem.Mol],
    fp_size: int = 2048,
    radius: int = 3,
    use_features: bool = True,
    as_array: bool = True,
) -> Optional[Union[np.ndarray, DataStructs.ExplicitBitVect]]:
    """Transform a molecule from smiles to morgan fingerprint.

    Args:
        mol (Chem.Mol or str): a molecule or a SMILES.
        fp_size (int, optional): Size of morgan fingerprint. Default to 2048.
        radius (int, optional): Radius of the morgan fingerprints. Default to 3.
        use_features: Whether to use atom features. Default to True.
        as_array: Whether to return a numpy array of an RDKit vec. Default to True.

    Returns:
        A fingerprint vector or None
    """

    if isinstance(mol, str):
        mol = dm.to_mol(mol)

    fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(
        mol,
        radius,
        nBits=fp_size,
        useFeatures=use_features,
    )

    if as_array:
        return fp_to_array(fp)

    return fp
