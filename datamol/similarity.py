from typing import List

from rdkit import Chem
from rdkit.DataManip.Metric import GetTanimotoDistMat

import numpy as np

import datamol as dm


def pdist(mols: List[Chem.Mol], **fp_args):
    """Compute the pairwise tanimoto distance between the fingerprints of all the
    molecules in the input set.

    Args:
        mols (list of Chem.Mol): list of molecules
        **fp_args (dict): list of args to pass to `to_fp()`.

    Returns:
        distmat, valid_idx: Distance matrix, and valid index that have passed the conversion
            to fingerprint.
    """

    # TODO(hadim): parallelize me.
    fps = [dm.to_fp(mol, as_array=False, **fp_args) for mol in mols]

    valid_idx, fps = zip(*[(i, fp) for i, fp in enumerate(fps) if fp is not None])
    fps = list(fps)

    dist = GetTanimotoDistMat(fps)
    dist_mat = np.zeros((len(fps), len(fps)))
    dist_mat[np.triu_indices_from(dist_mat, 1)] = dist
    dist_mat += dist_mat.T

    return dist_mat, np.array(valid_idx)
