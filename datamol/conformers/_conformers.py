from typing import Union
from typing import List
from typing import Dict
from typing import Any
from typing import Sequence
from typing import Optional

import copy

from loguru import logger

import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolTransforms
from rdkit.ML.Cluster import Butina

import datamol as dm


def generate(
    mol: Chem.rdchem.Mol,
    n_confs: int = None,
    rms_cutoff: Optional[float] = None,
    clear_existing: bool = True,
    align_conformers: bool = True,
    minimize_energy: bool = False,
    method: str = None,
    energy_iterations: int = 500,
    warning_not_converged: int = 10,
    random_seed: int = 19,
    add_hs: bool = True,
    verbose: bool = False,
) -> Chem.rdchem.Mol:
    """Compute conformers of a molecule.

    Example:

    ```python
    import datamol as dm
    smiles = "O=C(C)Oc1ccccc1C(=O)O"
    mol = dm.to_mol(smiles)
    mol = dm.conformers.generate(mol)

    # Get all conformers as a list
    conformers = mol.GetConformers()

    # Get the 3D atom positions of the first conformer
    positions = mol.GetConformer(0).GetPositions()

    # If minimization has been enabled (default to True)
    # you can access the computed energy.
    conf = mol.GetConformer(0)
    props = conf.GetPropsAsDict()
    print(props)
    # {'rdkit_uff_energy': 1.7649408317784008}
    ```

    Args:
        mol: a molecule
        n_confs: Number of conformers to generate. Depends on the
            number of rotatable bonds by default.
        rms_cutoff: The minimum RMS value in Angstrom at which two conformers
            are considered redundant and one is deleted. If None, all conformers
            are kept. This step is done after an eventual minimization step.
        clear_existing: Whether to overwrite existing conformers for the molecule.
        align_conformers: Wehther to align conformer.
        minimize_energy: Wether to minimize conformer's energies using UFF.
            Disable to generate conformers much faster.
        method: RDKit method to use for embedding. Choose among
            ["ETDG", "ETKDG", "ETKDGv2", "ETKDGv3"]. If None, "ETKDGv3" is used.
        energy_iterations: Maximum number of iterations during the energy minimization procedure.
            It corresponds to the `maxIters` argument in RDKit.
        warning_not_converged: Wether to log a warning when the number of not converged conformers
            during the minimization is higher than `warning_not_converged`. Only works when `verbose` is set to True. Disable with 0. Defaults to 10.
        random_seed: Set to None or -1 to disable.
        add_hs: Whether to add hydrogens to the mol before embedding. If set to True, the hydrogens
            are removed in the returned molecule. Warning: explicit hydrogens won't be conserved. It is strongly
            recommended to let the default value to True. The RDKit documentation says: "To get good 3D conformations,
            it’s almost always a good idea to add hydrogens to the molecule first."
        verbose: Wether to enable logs during the process.

    Returns:
        mol: the molecule with the conformers.
    """

    AVAILABLE_METHODS = ["ETDG", "ETKDG", "ETKDGv2", "ETKDGv3"]

    if method is None:
        method = "ETKDGv3"

    if method not in AVAILABLE_METHODS:
        raise ValueError(f"The method {method} is not supported. Use from {AVAILABLE_METHODS}")

    # Random seed
    if random_seed is None:
        random_seed = -1

    # Clone molecule
    mol = copy.deepcopy(mol)

    # Remove existing conformers
    if clear_existing:
        mol.RemoveAllConformers()

    # Add hydrogens
    if add_hs:
        mol = Chem.AddHs(mol)

    if not n_confs:
        # Set the number of conformers depends on
        # the number of rotatable bonds.
        rotatable_bonds = Descriptors.NumRotatableBonds(mol)
        if rotatable_bonds < 8:
            n_confs = 50
        elif rotatable_bonds < 12:
            n_confs = 200
        else:
            n_confs = 300

    # Embed conformers
    params = getattr(AllChem, method)()
    params.randomSeed = random_seed
    params.enforceChirality = True
    confs = AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, params=params)

    # Sometime embedding fails. Here we try again by disabling `enforceChirality`.
    if len(confs) == 0:
        if verbose:
            logger.warning(
                f"Conformers embedding failed for {dm.to_smiles(mol)}. Trying without enforcing chirality."
            )
        params = getattr(AllChem, method)()
        params.randomSeed = random_seed
        params.enforceChirality = False
        confs = AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, params=params)

    if len(confs) == 0:
        raise ValueError(f"Conformers embedding failed for {dm.to_smiles(mol)}")

    # Minimize energy
    if minimize_energy:

        # Minimize conformer's energy using UFF
        results = AllChem.UFFOptimizeMoleculeConfs(mol, maxIters=energy_iterations)
        energies = [energy for _, energy in results]

        # Some conformers might not have converged during minimization.
        not_converged = sum([not_converged for not_converged, _ in results if not_converged])
        if warning_not_converged != 0 and not_converged > warning_not_converged and verbose:
            logger.warning(
                f"{not_converged}/{len(results)} conformers have not converged for {dm.to_smiles(mol)}"
            )

        # Add the energy as a property to each conformers
        [
            conf.SetDoubleProp("rdkit_uff_energy", energy)
            for energy, conf in zip(energies, mol.GetConformers())
        ]

        # Now we reorder conformers according to their energies,
        # so the lowest energies conformers are first.
        mol_clone = copy.deepcopy(mol)
        ordered_conformers = [conf for _, conf in sorted(zip(energies, mol_clone.GetConformers()))]
        mol.RemoveAllConformers()
        [mol.AddConformer(conf, assignId=True) for conf in ordered_conformers]

    # Align conformers to each others
    if align_conformers:
        rdMolAlign.AlignMolConformers(mol)

    if rms_cutoff is not None:
        mol = cluster(
            mol,
            rms_cutoff=rms_cutoff,
            already_aligned=align_conformers,
            centroids=True,
        )  # type: ignore

    if add_hs:
        mol = Chem.RemoveHs(mol)

    return mol


def cluster(
    mol: Chem.rdchem.Mol,
    rms_cutoff: float = 1,
    already_aligned: bool = False,
    centroids: bool = True,
):
    """Cluster the conformers of a molecule according to an RMS threshold in Angstrom.

    Args:
        mol: a molecule
        rms_cutoff: The RMS cutoff in Angstrom.
        already_aligned: Whether or not the conformers are aligned. If False,
            they will be aligmned furing the RMS computation.
        centroids: If True, return one molecule with centroid conformers
            only. If False return a list of molecules per cluster with all
            the conformers of the cluster. Defaults to True.
    """

    # Clone molecule
    mol = copy.deepcopy(mol)

    # Compute RMS
    dmat = AllChem.GetConformerRMSMatrix(mol, prealigned=already_aligned)

    # Cluster
    conf_clusters = Butina.ClusterData(
        dmat,
        nPts=mol.GetNumConformers(),
        distThresh=rms_cutoff,
        isDistData=True,
        reordering=False,
    )

    return return_centroids(mol, conf_clusters, centroids=centroids)


def rmsd(mol: Chem.rdchem.Mol) -> np.ndarray:
    """Compute the RMSD between all the conformers of a molecule.

    Args:
        mol: a molecule
    """

    if mol.GetNumConformers() <= 1:
        raise ValueError(
            "The molecule has 0 or 1 conformer. You can generate conformers with `dm.conformers.generate(mol)`."
        )

    n_confs = mol.GetNumConformers()
    rmsds = []
    for i in range(n_confs):
        for j in range(n_confs):
            rmsd = rdMolAlign.AlignMol(prbMol=mol, refMol=mol, prbCid=i, refCid=j)
            rmsds.append(rmsd)
    return np.array(rmsds).reshape(n_confs, n_confs)


def return_centroids(
    mol: Chem.rdchem.Mol,
    conf_clusters: Sequence[Sequence[int]],
    centroids: bool = True,
) -> Union[List[Chem.rdchem.Mol], Chem.rdchem.Mol]:
    """Given a list of cluster indices, return one single molecule
    with only the centroid of the clusters of a list of molecules per cluster.

    Args:
        mol: a molecule.
        conf_clusters: list of cluster indices.
        centroids: If True, return one molecule with centroid conformers
            only. If False return a list of molecules per cluster with all
            the conformers of the cluster.
    """

    if centroids:
        # Collect centroid of each cluster (first element of the list)
        centroid_list = [indices[0] for indices in conf_clusters]

        # Keep only centroid conformers
        mol_clone = copy.deepcopy(mol)
        confs = [mol_clone.GetConformers()[i] for i in centroid_list]
        mol.RemoveAllConformers()
        [mol.AddConformer(conf, assignId=True) for conf in confs]
        return mol

    else:
        # Create a new molecule for each cluster and add conformers to it.
        mols = []
        for cluster in conf_clusters:
            m = copy.deepcopy(mol)
            m.RemoveAllConformers()
            [m.AddConformer(mol.GetConformer(c), assignId=True) for c in cluster]
            mols.append(m)
        return mols


def translate(mol: Chem.rdchem.Mol, new_centroid: Union[np.ndarray, List[int]], conf_id: int = -1):
    """Move a given conformer of a molecule to a new position. The transformation is performed
    in place.

    Args:
        mol: the molecule.
        new_centroid: the new position to move to of shape [x, y, z]
        conf_id: id of the conformer.
    """

    # Get conformer
    conf = mol.GetConformer(conf_id)

    # Compute the vector for translation
    mol_center = rdMolTransforms.ComputeCentroid(conf)
    mol_center = np.array([mol_center.x, mol_center.y, mol_center.z])

    # Make the transformation matrix
    T = np.eye(4)
    T[:3, 3] = new_centroid - mol_center

    # Transform
    rdMolTransforms.TransformConformer(conf, T)
