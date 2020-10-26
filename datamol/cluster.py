from typing import List
from typing import Callable
from typing import Optional

import operator
import functools

import numpy as np

from rdkit import Chem
from rdkit.Chem import DataStructs
from rdkit.ML.Cluster import Butina
from rdkit.SimDivFilters.rdSimDivPickers import MaxMinPicker

import datamol as dm


def cluster_mols(
    mols: List[Chem.Mol],
    cutoff: float = 0.2,
    feature_fn: Callable = None,
    n_jobs: Optional[int] = 1,
):
    """Cluster a set of molecules using the butina clustering algorithm and a given threshold.

    Args:
        mols (List[Chem.Mol]): a list of molecules.
        cutoff (float, optional): Cuttoff for the clustering. Default to 0.2.
        feature_fn (callable, optional): A feature function that takes a Chem.Mol object
            and return molecular features. By default, the `dm.to_fp()` is used.
            Default to None.
        n_jobs (int): Number of jobs for parallelization. Let to 1 for no
            parallelization. Set to None to use all available cores.
    """

    if feature_fn is None:
        feature_fn = functools.partial(dm.to_fp, as_array=False)

    features = dm.parallelized(feature_fn, mols, n_jobs=n_jobs)

    dists = []
    n_mols = len(mols)

    for i in range(1, n_mols):
        dist = DataStructs.BulkTanimotoSimilarity(features[i], features[:i], returnDistance=True)
        dists.extend([x for x in dist])

    # now cluster the data
    cluster_indices = Butina.ClusterData(dists, n_mols, cutoff, isDistData=True)
    cluster_mols = [operator.itemgetter(*cluster)(mols) for cluster in cluster_indices]

    # Make single mol cluster a list
    cluster_mols = [[c] if isinstance(c, Chem.Mol) else c for c in cluster_mols]

    return cluster_indices, cluster_mols


def pick_diverse(
    mols: List[Chem.Mol],
    npick: int,
    initial_picks: List[Chem.Mol] = None,
    feature_fn: Callable = None,
    dist_fn: Callable = None,
    seed: int = 42,
    n_jobs: Optional[int] = 1,
):
    r"""Pick a set of diverse molecules based on they fingerprint.

    Args:
        npick (int): Number of element to pick from mols.
        initial_picks (list, optional): Starting list of molecules that should be in the
            set of picked molecules. Default to None.
        feature_fn (callable, optional): A feature function that takes a Chem.Mol object
            and return molecular features. By default, the `dm.to_fp()` is used.
            Default to None.
        dist_fn (callable, optional): A function that takes two indexes (i,j) and return the
            distance between them. You might use partial to set the fingerprints as input.
            By default, the Tanimoto similarity will be used. Default to None.
        n_jobs (int): Number of jobs for parallelization. Let to 1 for no
            parallelization. Set to None to use all available cores.

    Returns:
        picked_inds: index of the molecule that have been picked
        mols: molecules that have been picked
    """

    if feature_fn is None:
        feature_fn = functools.partial(dm.to_fp, as_array=False)

    features = dm.parallelized(feature_fn, mols, n_jobs=n_jobs)

    def distij(i, j, features=features):
        return 1.0 - DataStructs.TanimotoSimilarity(features[i], features[j])

    if dist_fn is None:
        dist_fn = distij

    picker = MaxMinPicker()
    initial_picks = [] if initial_picks is None else initial_picks
    picked_inds = picker.LazyPick(dist_fn, len(mols), npick, firstPicks=initial_picks, seed=seed)
    picked_inds = np.array(picked_inds)
    picked_mols = [mols[x] for x in picked_inds]

    return picked_inds, picked_mols
