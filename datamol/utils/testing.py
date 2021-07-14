from typing import List
from typing import Optional
from typing import Union

import functools

import numpy as np
from scipy.spatial import distance

from rdkit import Chem
from rdkit.DataManip.Metric import GetTanimotoDistMat  # type: ignore
from rdkit.DataStructs.cDataStructs import TanimotoSimilarity

import datamol as dm


def pdist_rdkit(
    mols: List[Union[str, Chem.rdchem.Mol]],
    n_jobs: Optional[int] = 1,
    squareform: bool = True,
    **fp_args,
) -> np.ndarray:
    """Equivalent to `dm.similarity.pdist` but uses the RDKit API.

    Important:
        This function is only used for testing and shoult not be used in production.
    """

    fps = dm.parallelized(
        functools.partial(dm.to_fp, as_array=False, **fp_args),
        mols,
        n_jobs=n_jobs,
    )

    fps = list(fps)  # type: ignore

    dist = GetTanimotoDistMat(fps)

    # Put in squareform: `scipy.spatial.distance.squareform` is incompatible with RDKit returned vector.
    dist_mat = np.zeros((len(fps), len(fps)))
    dist_mat[np.tril_indices_from(dist_mat, -1)] = dist
    dist_mat += dist_mat.T

    if not squareform:
        dist_mat = distance.squareform(dist_mat, force="tovector")

    return dist_mat


def cdist_rdkit(
    mols1: List[Union[str, Chem.rdchem.Mol]],
    mols2: List[Union[str, Chem.rdchem.Mol]],
    n_jobs: Optional[int] = 1,
    **fp_args,
) -> np.ndarray:
    """Equivalent to `dm.similarity.cdist` but uses the RDKit API.

    Important:
        This function is only used for testing and shoult not be used in production.
    """

    fps1 = dm.parallelized(
        functools.partial(dm.to_fp, as_array=False, **fp_args),
        mols1,
        n_jobs=n_jobs,
    )

    fps2 = dm.parallelized(
        functools.partial(dm.to_fp, as_array=False, **fp_args),
        mols2,
        n_jobs=n_jobs,
    )

    fps1 = list(fps1)  # type: ignore
    fps2 = list(fps2)  # type: ignore

    dist_mat = np.zeros((len(fps1), len(fps2)))
    for i in range(len(fps1)):
        for j in range(len(fps2)):
            d = 1 - TanimotoSimilarity(fps1[i], fps2[j])
            dist_mat[i, j] = d

    return dist_mat
