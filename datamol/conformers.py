import copy

from loguru import logger

import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdFreeSASA
from rdkit.Chem import TorsionFingerprints
from rdkit.ML.Cluster import Butina


def generate(
    mol: Chem.Mol,
    n_confs: int = None,
    random_seed: int = 19,
    method: str = "ETKDGv3",
    do_clustering: bool = True,
    cluster_mode: str = "RMS",
    distance_threshold: float = 1,
    max_iterations: int = 500,
    warning_not_converged: int = 10,
    return_centroids: bool = True,
):
    """Compute conformers of a molecule.

    Args:
        mol (Chem.Mol): a molecule
        n_confs (int, optional): Number of conformers to generate. Depends on the
            number of rotatable bonds by default. Defaults to None.
        random_seed (int, optional): Set to None or -1 to disable. Defaults to 19.
        do_clustering (bool, optional): Whether to do clustering or not.
        cluster_mode (str, optional): Clustering mode. Use "RMS" or "TFD". Defaults to "RMS".
        distance_threshold (float, optional): Threshold for the clustering. Defaults to 1.
        max_iterations (int, optional): Maximum number of iterations during UFF minimization.
            Defaults to 500.
        warning_not_converged (int, optional): Display a warning if n conformers did not
            converged during minimization. Defaults to 10.
        return_centroids (bool, optional): If True, return one molecule with centroid conformers
            only sorted by lowest energy. If False return a list of molecules per cluster with all
            the cluster conformers. Defaults to True.
    Returns:
        [type]: [description]
    """

    AVAILABLE_METHODS = ["ETDG", "ETKDG", "ETKDGv2", "ETKDGv3"]

    # Random seed
    if random_seed is None:
        random_seed = -1

    # Clone molecule
    mol = copy.deepcopy(mol)

    # It's probably best to compute conformers using the hydrogen
    # atoms. Even if rdkit does not always use those.
    mol = Chem.AddHs(mol)

    if not n_confs:
        # Set the number of conformers depends on
        # the number of rotatable bonds.
        rotatable_bonds = Chem.Descriptors.NumRotatableBonds(mol)
        if rotatable_bonds < 8:
            n_confs = 50
        elif rotatable_bonds < 12:
            n_confs = 200
        else:
            n_confs = 300

    # Embed conformers
    params = AllChem.ETKDGv2()
    params.randomSeed = random_seed
    params.enforceChirality = True
    confs = AllChem.EmbedMultipleConfs(clone, numConfs=n_confs, params=params)

    # Sometime embedding fails. Here we try again by disabling `enforceChirality`.
    if len(confs) == 0:
        logger.warning(
            f"Conformers embedding failed for {Chem.MolToSmiles(mol)}. Trying without enforcing chirality."
        )
        params = AllChem.ETKDGv2()
        params.randomSeed = random_seed
        params.enforceChirality = False
        confs = AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, params=params)

    assert len(confs) > 0, f"Conformers embedding failed for {Chem.MolToSmiles(mol)}"

    # Align conformers to each others
    Chem.rdMolAlign.AlignMolConformers(mol)

    # Minimize energy
    results = AllChem.UFFOptimizeMoleculeConfs(mol, maxIters=max_iterations)

    energies = [energy for _, energy in results]
    not_converged = sum([not_converged for not_converged, _ in results if not_converged])
    if not_converged > warning_not_converged:
        logger.warning(
            f"{not_converged}/{len(results)} conformers have not converged for {Chem.MolToSmiles(mol)}"
        )

    return mol


def cluster(
    mol: Chem.Mol,
    n_confs: int = None,
    random_seed: int = 19,
    method: str = "ETKDGv3",
    do_clustering: bool = True,
    cluster_mode: str = "RMS",
    distance_threshold: float = 1,
    max_iterations: int = 500,
    warning_not_converged: int = 10,
    return_centroids: bool = True,
):
    """Compute conformers of a molecule.

    Args:
        mol (Chem.Mol): a molecule
        n_confs (int, optional): Number of conformers to generate. Depends on the
            number of rotatable bonds by default. Defaults to None.
        random_seed (int, optional): Set to None to disable. Defaults to 19.
        do_clustering (bool, optional): Whether to do clustering or not.
        cluster_mode (str, optional): Clustering mode. Use "RMS" or "TFD". Defaults to "RMS".
        distance_threshold (float, optional): Threshold for the clustering. Defaults to 1.
        max_iterations (int, optional): Maximum number of iterations during UFF minimization.
            Defaults to 500.
        warning_not_converged (int, optional): Display a warning if n conformers did not
            converged during minimization. Defaults to 10.
        return_centroids (bool, optional): If True, return one molecule with centroid conformers
            only sorted by lowest energy. If False return a list of molecules per cluster with all
            the cluster conformers. Defaults to True.
    Returns:
        [type]: [description]
    """

    if not do_clustering:
        # Sort conformers by energies
        indices = np.argsort(energies)
        confs = [mol.GetConformers()[i] for i in indices]
        mol.RemoveAllConformers()
        [mol.AddConformer(conf, assignId=True) for conf in confs]
        return mol

    # Cluster conformers
    if cluster_mode == "TFD":
        dmat = TorsionFingerprints.GetTFDMatrix(
            mol, useWeights=False, maxDev="equal", symmRadius=2, ignoreColinearBonds=True
        )
    elif cluster_mode == "RMS":
        dmat = AllChem.GetConformerRMSMatrix(mol, prealigned=False)
    else:
        raise ValueError(f"{cluster_mode} wrong mode")

    conf_clusters = Butina.ClusterData(
        dmat, nPts=len(confs), distThresh=distance_threshold, isDistData=True, reordering=False
    )

    # Collect centroid of each cluster (first element of the list)
    centroids = [indices[0] for indices in conf_clusters]

    if return_centroids:
        # Get energy for centroids
        centroid_energies = [energies[i] for i in centroids]

        # Sort centroids by lowest energy
        centroid_indices = np.argsort(centroid_energies)

        # Keep only centroid conformers
        confs = [mol.GetConformers()[centroids[i]] for i in centroid_indices]
        mol.RemoveAllConformers()
        [mol.AddConformer(conf, assignId=True) for conf in confs]
        return Chem.RemoveHs(mol)

    else:
        # Create a new molecule for each cluster and add conformers to it.
        mols = []
        for cluster in conf_clusters:
            m = copy.deepcopy(mol)
            m.RemoveAllConformers()
            [m.AddConformer(mol.GetConformer(c), assignId=True) for c in cluster]
            mols.append(Chem.RemoveHs(m))
        return mols

def sasa(mol, conf_id=0):
    """Compute Solvent Accessible Surface Area.

    TODO(hadim): handle when no conformers are present.
    """
    mol = Chem.AddHs(mol)

    # Get Van der Waals radii (angstrom)
    ptable = Chem.GetPeriodicTable()
    radii = [ptable.GetRvdw(atom.GetAtomicNum()) for atom in mol.GetAtoms()]

    # Compute solvent accessible surface area
    sasa = rdFreeSASA.CalcSASA(mol, radii, confIdx=conf_id)

    if np.isnan(sasa):
        return None

    return sasa


def rmsd(mol):
    """TODO(hadim): handle when no conformers are present."""
    n_confs = mol.GetNumConformers()
    rmsds = []
    for i in range(n_confs):
        for j in range(n_confs):
            rmsd = Chem.rdMolAlign.AlignMol(prbMol=mol, refMol=mol, prbCid=i, refCid=j)
            rmsds.append(rmsd)
    return np.array(rmsds).reshape(n_confs, n_confs)
