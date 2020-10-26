from typing import List
from typing import Optional

import functools

from rdkit import Chem
from rdkit.DataManip.Metric import GetTanimotoDistMat

import numpy as np

import datamol as dm


def pdist(mols: List[Chem.Mol], n_jobs: Optional[int] = 1, **fp_args):
    """Compute the pairwise tanimoto distance between the fingerprints of all the
    molecules in the input set.

    Args:
        mols (list of Chem.Mol): list of molecules
        n_jobs (int): Number of jobs for parallelization. Let to 1 for no
            parallelization. Set to None to use all available cores.
        **fp_args (dict): list of args to pass to `to_fp()`.

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
