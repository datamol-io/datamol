from typing import List
from typing import Optional
from typing import Union
from typing import Any

import functools

import numpy as np
from sklearn.metrics import pairwise_distances_chunked
from scipy.spatial import distance

import datamol as dm


def pdist(
    mols: List[Union[str, dm.Mol]],
    n_jobs: Optional[int] = 1,
    squareform: bool = True,
    **fp_args: Any,
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
    mols1: List[Union[str, dm.Mol]],
    mols2: List[Union[str, dm.Mol]],
    n_jobs: Optional[int] = 1,
    distances_chunk: bool = False,
    distances_chunk_memory: int = 1024,
    distances_n_jobs: int = -1,
    **fp_args: Any,
) -> np.ndarray:
    """Compute the tanimoto distance between the fingerprints of each pair of
    molecules of the two collections of inputs.

    Args:
        mols1: list of molecules.
        mols2: list of molecules.
        n_jobs: Number of jobs for fingerprint computation. Let to 1 for no
            parallelization. Set to None or -1 to use all available cores.
        distances_chunk: Whether to use chunked computation.
        distances_chunk_memory: Memory size in MB to use for chunked computation.
        distances_n_jobs: Number of jobs for parallelization.
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

    fps1 = np.array(fps1).astype(bool)
    fps2 = np.array(fps2).astype(bool)

    if distances_chunk:
        distances = pairwise_distances_chunked(
            fps1,
            fps2,
            metric="jaccard",
            n_jobs=distances_n_jobs,
            working_memory=distances_chunk_memory,
        )
        distances = [i for i in distances]
        distances = np.vstack(distances)
    else:
        distances = distance.cdist(fps1, fps2, metric="jaccard")

    return distances
