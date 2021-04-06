from typing import List
from typing import Callable
from typing import Optional

import operator
import functools
from collections import defaultdict as ddict

from loguru import logger

import numpy as np
from scipy.spatial import distance

from rdkit import Chem
from rdkit.Chem import DataStructs
from rdkit.ML.Cluster import Butina

from rdkit.SimDivFilters.rdSimDivPickers import ClusterMethod
from rdkit.SimDivFilters.rdSimDivPickers import HierarchicalClusterPicker
from rdkit.SimDivFilters.rdSimDivPickers import MaxMinPicker
from rdkit.SimDivFilters.rdSimDivPickers import LeaderPicker

import datamol as dm


def cluster_mols(
    mols: List[Chem.rdchem.Mol],
    cutoff: float = 0.2,
    feature_fn: Callable = None,
    n_jobs: Optional[int] = 1,
):
    """Cluster a set of molecules using the butina clustering algorithm and a given threshold.

    Args:
        mols: a list of molecules.
        cutoff: Cuttoff for the clustering. Default to 0.2.
        feature_fn: A feature function that takes a Chem.rdchem.Mol object
            and return molecular features. By default, the `dm.to_fp()` is used.
            Default to None.
        n_jobs: Number of jobs for parallelization. Let to 1 for no
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
    cluster_mols = [[c] if isinstance(c, Chem.rdchem.Mol) else c for c in cluster_mols]

    return cluster_indices, cluster_mols


def pick_diverse(
    mols: List[Chem.rdchem.Mol],
    npick: int,
    initial_picks: List[int] = None,
    feature_fn: Callable = None,
    dist_fn: Callable = None,
    seed: int = 42,
    n_jobs: Optional[int] = 1,
):
    r"""Pick a set of diverse molecules based on they fingerprint.

    Args:
        mols: a list of molecules.
        npick: Number of element to pick from mols, including the preselection.
        initial_picks: Starting list of index for molecules that should be in the
            set of picked molecules. Default to None.
        feature_fn: A feature function that takes a Chem.rdchem.Mol object
            and return molecular features. By default, the `dm.to_fp()` is used.
            Default to None.
        dist_fn: A function that takes two indexes (i,j) and return the
            distance between them. You might use partial to set the fingerprints as input.
            By default, the Tanimoto similarity will be used. Default to None.
        seed: seed for reproducibility
        n_jobs: Number of jobs for parallelization. Let to 1 for no
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


def pick_centroids(
    mols: List[Chem.rdchem.Mol],
    npick: int = 0,
    initial_picks: List[int] = None,
    threshold: float = 0.5,
    feature_fn: Callable = None,
    dist_fn: Callable = None,
    seed: int = 42,
    method: str = "sphere",
    n_jobs: Optional[int] = 1,
):
    r"""Pick a set of `npick` centroids from a list of molecules.

    Args:
        mols: a list of molecules.
        npick: Number of element to pick from mols, including the preselection.
        threshold: Minimum distance between centroids for `maxmin` and sphere exclusion (`sphere`) methods.
        initial_picks: Starting list of index for molecules that should be in the
            set of picked molecules. Default to None.
        feature_fn (callable, optional): A feature function that takes a Chem.rdchem.Mol object
            and return molecular features. By default, the `dm.to_fp()` is used.
            Default to None.
        dist_fn: A function that takes two indexes (i,j) and return the
            distance between them. You might use partial to set the fingerprints as input.
            By default, the Tanimoto similarity will be used. Default to None.
        seed: seed for reproducibility
        method: Picking method to use. One of  `sphere`, `maxmin` or any
            supported rdkit hierarchical clustering method such as `centroid`, `clink`, `upgma`
        n_jobs: Number of jobs for parallelization. Let to 1 for no
            parallelization. Set to None to use all available cores.

    Returns:
        picked_inds: index of the molecule that have been selected as centroids
        mols: molecules that have been picked
    """

    n_mols = len(mols)
    if feature_fn is None:
        feature_fn = functools.partial(dm.to_fp, as_array=False)

    features = dm.parallelized(feature_fn, mols, n_jobs=n_jobs)

    def distij(i, j, features=features):
        return 1.0 - DataStructs.TanimotoSimilarity(features[i], features[j])

    if dist_fn is None:
        dist_fn = distij

    initial_picks = [] if initial_picks is None else initial_picks

    if method == "maxmin":
        picker = MaxMinPicker()
        picked_inds, _ = picker.LazyPickWithThreshold(
            dist_fn,
            n_mols,
            pickSize=npick,
            threshold=threshold,
            firstPicks=initial_picks,
            seed=seed,
        )

    elif method == "sphere":
        picker = LeaderPicker()
        picked_inds = picker.LazyPick(
            dist_fn, n_mols, threshold=threshold, pickSize=npick, firstPicks=initial_picks
        )

    elif method.upper() in ClusterMethod.names.keys() and npick:
        if initial_picks:
            logger.warning(
                "Initial picks is not supported by hierarchical clustering. You pick has been discarded."
            )

        dist_mat = dm.parallelized(
            distij, list(zip(*np.tril_indices(len(mols), k=-1))), arg_type="args"
        )
        dist_mat = np.asarray(dist_mat)
        picker = HierarchicalClusterPicker(ClusterMethod.names[method.upper()])
        picked_inds = picker.Pick(dist_mat, n_mols, npick)
    else:
        raise ValueError(f"Picking method {method} with {npick} elements to pick is not supported.")
    picked_inds = np.array(picked_inds)
    picked_mols = [mols[x] for x in picked_inds]

    return picked_inds, picked_mols


def assign_to_centroids(
    mols: List[Chem.rdchem.Mol],
    centroids: List[Chem.rdchem.Mol],
    feature_fn: Callable = None,
    dist_fn: Callable = None,
    n_jobs: Optional[int] = 1,
):
    r"""Assign molecules to centroids. Each molecule will be assigned to the closest centroid.

    Args:
        mols: a list of molecules to assign to centroids
        centroids: list of molecules to use as centroid
        feature_fn: A feature function that takes a Chem.rdchem.Mol object
            and return molecular features. By default, the `dm.to_fp()` is used.
            Default to None.
        dist_fn: A function that takes two indexes (i,j) and return the
            distance between them. You might use partial to set the fingerprints as input.
            By default, the Tanimoto similarity will be used. Default to None.
        n_jobs: Number of jobs for parallelization. Let to 1 for no
            parallelization. Set to None to use all available cores.

    Returns:
        clusters_map: dict of index mapping each centroid index to the molecule index in the cluster
        clusters_list: list of all molecules in each cluster. The cluster index follows the index of the centroid.
            Note that the centroid molecule is not added to the cluster.
    """

    if feature_fn is None:
        feature_fn = functools.partial(dm.to_fp, as_array=False)

    all_mols = [x for x in mols] + [c for c in centroids]
    features = dm.parallelized(feature_fn, all_mols, n_jobs=n_jobs)

    def distij(i, j, features=features):
        return 1.0 - DataStructs.TanimotoSimilarity(features[int(i)], features[int(j)])

    if dist_fn is None:
        dist_fn = distij

    clusters_map = ddict(list)
    clusters_list = [[] for _ in centroids]
    query_inds = np.expand_dims(np.arange(len(mols), dtype=int), axis=1)
    centroid_inds = np.expand_dims(np.arange(len(centroids), dtype=int), axis=1) + len(mols)
    dist_mat = distance.cdist(query_inds, centroid_inds, metric=distij)
    closest = np.argmin(dist_mat, axis=1)
    for ind, cluster_ind in enumerate(closest):  # type: ignore
        clusters_map[cluster_ind].append(ind)
        clusters_list[cluster_ind].append(mols[ind])
    return clusters_map, clusters_list
