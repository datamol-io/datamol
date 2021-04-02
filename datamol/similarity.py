from typing import List
from typing import Optional
from typing import Tuple

import functools

from rdkit import Chem
from rdkit.DataManip.Metric import GetTanimotoDistMat

import numpy as np
from scipy.spatial import distance

import datamol as dm


def pdist(
    mols: List[Chem.rdchem.Mol], n_jobs: Optional[int] = 1, **fp_args
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute the pairwise tanimoto distance between the fingerprints of all the
    molecules in the input set.

    Args:
        mols: list of molecules
        n_jobs: Number of jobs for parallelization. Let to 1 for no
            parallelization. Set to None to use all available cores.
        **fp_args: list of args to pass to `to_fp()`.

    Returns:
        distmat, valid_idx: Distance matrix, and valid index that have passed the conversion
            to fingerprint.
    """

    fps = dm.parallelized(
        functools.partial(dm.to_fp, as_array=False, **fp_args),
        mols,
        n_jobs=n_jobs,
    )

    valid_idx, fps = zip(*[(i, fp) for i, fp in enumerate(fps) if fp is not None])
    fps = list(fps)

    dist = GetTanimotoDistMat(fps)
    dist_mat = np.zeros((len(fps), len(fps)))
    dist_mat[np.triu_indices_from(dist_mat, 1)] = dist
    dist_mat += dist_mat.T

    return dist_mat, np.array(valid_idx)


def cdist(
    mols1: List[Chem.rdchem.Mol],
    mols2: List[Chem.rdchem.Mol],
    n_jobs: Optional[int] = 1,
    **fp_args,
) -> np.ndarray:
    """Compute the pairwise tanimoto distance between the fingerprints of
    each pair of molecules of the two collections of inputs.

    Args:
        mols1: list of molecules.
        mols2: list of molecules.
        n_jobs: Number of jobs for parallelization. Let to 1 for no
            parallelization. Set to None to use all available cores.
        **fp_args: list of args to pass to `to_fp()`.

    Returns:
        distmat
    """

    fps1 = dm.parallelized(
        functools.partial(dm.to_fp, as_array=True, **fp_args),
        mols1,
        n_jobs=n_jobs,
    )

    fps2 = dm.parallelized(
        functools.partial(dm.to_fp, as_array=True, **fp_args),
        mols2,
        n_jobs=n_jobs,
    )

    fps1 = np.array(fps1)
    fps2 = np.array(fps2)

    dist_mat = distance.cdist(fps1, fps2, metric="jaccard")

    return dist_mat
