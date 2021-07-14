from typing import List
from typing import Optional
from typing import Union

import functools

from rdkit import Chem

import numpy as np
from scipy.spatial import distance

import datamol as dm


def pdist(
    mols: List[Union[str, Chem.rdchem.Mol]],
    n_jobs: Optional[int] = 1,
    squareform: bool = True,
    **fp_args,
) -> np.ndarray:
    """Compute the pairwise tanimoto distance between the fingerprints of all the
    molecules in the input set.

    Args:
        mols: list of molecules
        n_jobs: Number of jobs for parallelization. Let to 1 for no
            parallelization. Set to None to use all available cores.
        squareform: Whether to return in square form (matrix) or in a condensed
            form (1D vector).
        **fp_args: list of args to pass to `to_fp()`.

    Returns:
        dist_mat
    """

    fps = dm.parallelized(
        functools.partial(dm.to_fp, as_array=True, **fp_args),
        mols,
        n_jobs=n_jobs,
    )

    fps = np.array(fps)

    dist_mat = distance.pdist(fps, metric="jaccard")

    if squareform:
        dist_mat = distance.squareform(dist_mat, force="tomatrix")

    return dist_mat


def cdist(
    mols1: List[Union[str, Chem.rdchem.Mol]],
    mols2: List[Union[str, Chem.rdchem.Mol]],
    n_jobs: Optional[int] = 1,
    **fp_args,
) -> np.ndarray:
    """Compute the tanimoto distance between the fingerprints of each pair of
    molecules of the two collections of inputs.

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
