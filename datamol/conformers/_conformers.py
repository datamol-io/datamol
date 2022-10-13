from typing import Union
from typing import List
from typing import Sequence
from typing import Optional
from typing import Tuple

import copy

from loguru import logger

import numpy as np

from rdkit.Chem import AllChem
from rdkit.Chem import rdDistGeom
from rdkit.Chem import rdMolAlign
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdMolTransforms
from rdkit.Chem import rdForceFieldHelpers
from rdkit.ML.Cluster import Butina

from .. import Mol
from .. import mol as dm_mol
from .. import convert
from .. import descriptors


def generate(
    mol: Mol,
    n_confs: Optional[int] = None,
    use_random_coords: bool = True,
    enforce_chirality: bool = True,
    num_threads: int = 1,
    rms_cutoff: Optional[float] = None,
    clear_existing: bool = True,
    align_conformers: bool = True,
    minimize_energy: bool = False,
    sort_by_energy: bool = True,
    method: Optional[str] = None,
    forcefield: str = "UFF",
    ewindow: float = np.inf,
    eratio: float = np.inf,
    energy_iterations: int = 200,
    warning_not_converged: int = 0,
    random_seed: int = 19,
    add_hs: bool = True,
    ignore_failure: bool = False,
    embed_params: Optional[dict] = None,
    verbose: bool = False,
) -> Mol:
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
    conf = mol.GetConformer(1)
    props = conf.GetPropsAsDict()
    print(props)
    # {'rdkit_UFF_energy': 35.64074017773132,'rdkit_UFF_delta_energy': 0.24682258222552633}
    ```

    Args:
        mol: a molecule
        n_confs: Number of conformers to generate. Depends on the number of rotatable bonds
            by default: 50 for <8, 200 for <12 and 300 for >12.
        use_random_coords: Start the embedding from random coordinates instead of using eigenvalues
            of the distance matrix.
        enforce_chirality: Enforce correct chirilaty if chiral centers are present.
        num_threads: Number of threads to use when embedding multiple conformations.
        rms_cutoff: The minimum RMS value in Angstrom at which two conformers
            are considered redundant and one is deleted. If None, all conformers
            are kept. This step is done after an eventual minimization step.
        clear_existing: Whether to overwrite existing conformers for the molecule.
        align_conformers: Whether to align the conformers.
        minimize_energy: Whether to minimize conformer's energies using MMFF94s.
            Disable to generate conformers much faster.
        sort_by_energy: Sort conformers by energy when minimizing is turned to False.
        method: RDKit method to use for embedding. Choose among
            ["ETDG", "ETKDG", "ETKDGv2", "ETKDGv3"]. If None, "ETKDGv3" is used.
        forcefield: molecular forcefield to use, one of ['UFF','MMFF94s','MMFF94s_noEstat']
        ewindow: maximum energy above minimum energy conformer to output
        eratio: max delta-energy divided by rotatable bonds for conformers
        energy_iterations: Maximum number of iterations during the energy minimization procedure.
            It corresponds to the `maxIters` argument in RDKit.
        warning_not_converged: Wether to log a warning when the number of not converged conformers
            during the minimization is higher than `warning_not_converged`. Only works when `verbose` is set to True.
            Disable with 0. Defaults to 10.
        random_seed: Set to None or -1 to disable.
        add_hs: Whether to add hydrogens to the mol before embedding. If set to True, the hydrogens
            are removed in the returned molecule. Warning: explicit hydrogens won't be conserved. It is strongly
            recommended to let the default value to True. The RDKit documentation says: "To get good 3D conformations,
            it's almost always a good idea to add hydrogens to the molecule first."
        ignore_failure: It set to True, this will avoid raising an error when the embedding fails and return None instead.
        embed_params: Allows the user to specify arbitrary embedding parameters for the conformers. This will override any
            other default settings. See https://www.rdkit.org/docs/source/rdkit.Chem.rdDistGeom.html#rdkit.Chem.rdDistGeom.EmbedParameters
            for more details.
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
        mol = dm_mol.add_hs(mol)

    rotatable_bonds = descriptors.n_rotatable_bonds(mol)
    if not n_confs:
        # Set the number of conformers depends on
        # the number of rotatable bonds.
        if rotatable_bonds < 8:
            n_confs = 50
        elif rotatable_bonds < 12:
            n_confs = 200
        else:
            n_confs = 300

    # Setup the parameters for the embedding
    params = getattr(rdDistGeom, method)()
    params.randomSeed = random_seed
    params.enforceChirality = enforce_chirality
    params.useRandomCoords = use_random_coords
    params.numThreads = num_threads

    if embed_params is not None:
        for k, v in embed_params.items():
            setattr(params, k, v)

    # Embed conformers
    confs = rdDistGeom.EmbedMultipleConfs(mol, numConfs=n_confs, params=params)

    if len(confs) == 0:
        if ignore_failure:
            if verbose:
                logger.warning(
                    f"Conformers embedding failed for {convert.to_smiles(mol)}. Returning None because ignore_failure is set."
                )
            return None
        raise ValueError(f"Conformers embedding failed for {convert.to_smiles(mol)}")

    energies = None

    # Minimize energy
    if minimize_energy:

        # Minimize conformer's energy using MMFF
        ff = _get_ff(mol, forcefield)
        results = rdForceFieldHelpers.OptimizeMoleculeConfs(
            mol, ff, maxIters=energy_iterations, numThreads=num_threads
        )
        energies = [energy for _, energy in results]

        # Some conformers might not have converged during minimization.
        not_converged = sum([not_converged for not_converged, _ in results if not_converged])
        if warning_not_converged != 0 and not_converged > warning_not_converged and verbose:
            logger.warning(
                f"{not_converged}/{len(results)} conformers have not converged for {convert.to_smiles(mol)}"
            )

    elif sort_by_energy:
        energies = []
        for conf in mol.GetConformers():
            ff = _get_ff(mol, forcefield, conf_id=conf.GetId())
            energies.append(ff.CalcEnergy())
        energies = np.array(energies)

    if energies is not None:
        minE = np.min(energies)
        # Add the energy as a property to each conformers
        [
            (
                conf.SetDoubleProp(f"rdkit_{forcefield}_energy", energy),
                conf.SetDoubleProp(f"rdkit_{forcefield}_delta_energy", energy - minE),
            )
            for energy, conf in zip(energies, mol.GetConformers())
        ]

        # Now we reorder conformers according to their energies,
        # so the lowest energies conformers are first.  eliminate conformers that exceed ewindow, eratio
        mol_clone = copy.deepcopy(mol)
        ordered_conformers = [
            conf
            for E, conf in sorted(zip(energies, mol_clone.GetConformers()), key=lambda x: x[0])
            if E - minE <= ewindow and (E - minE) / rotatable_bonds <= eratio
        ]
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
        mol = dm_mol.remove_hs(mol)

    return mol


def _get_ff(mol: Mol, forcefield: str, conf_id: int = -1):
    """Gets molecular forcefield for input mol according to name
    Args:
        mol: input molecule
        forcefield: forcefield name. One of "UFF", "MMFF94s", "MMFF94s_noEstat"]
        conf_id: conformer id. -1 is used by default
    """
    assert forcefield in [
        "UFF",
        "MMFF94s",
        "MMFF94s_noEstat",
    ], f"Forcefield {forcefield} is not supported"
    if forcefield == "UFF":
        return rdForceFieldHelpers.UFFGetMoleculeForceField(mol, confId=conf_id)

    mp = rdForceFieldHelpers.MMFFGetMoleculeProperties(mol, "MMFF94s")
    if forcefield == "MMFF94s_noEstat":
        mp.SetMMFFEleTerm(False)
    return rdForceFieldHelpers.MMFFGetMoleculeForceField(mol, mp, confId=conf_id)


def cluster(
    mol: Mol,
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


def rmsd(mol: Mol) -> np.ndarray:
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
    mol: Mol,
    conf_clusters: Sequence[Sequence[int]],
    centroids: bool = True,
) -> Union[List[Mol], Mol]:
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


def translate(mol: Mol, new_centroid: Union[np.ndarray, List[int]], conf_id: int = -1):
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


def align_conformers(
    mols: List[Mol],
    ref_id: int = 0,
    copy: bool = True,
    conformer_id: int = -1,
    backend: str = "crippenO3A",
) -> Tuple[list, list]:
    """Align a list of molecules to a reference molecule.

    Note that using the `O3A` backend, hydrogens will be added at the beginning of the procedure
    and removed at the end of the procedure.

    Args:
        mols: List of molecules to align. All the molecules must have a conformer.
        ref_id: Index of the reference molecule. By default, the first molecule in the list
            will be used as reference.
        copy: Whether to copy the molecules before performing the alignement.
        conformer_id: Conformer id to use.
        backend: Backend to use to compute the alignment from `crippenO3A`, `O3A`.

    Returns:
        mols: The aligned molecules.
        scores: The score of the alignement.
    """

    allowed_backends = ["crippenO3A", "O3A"]
    if backend not in allowed_backends:
        raise ValueError(
            f"The backend '{backend}' is not supported. Choose from: {allowed_backends}"
        )

    # Check all input molecules has a conformer
    if not all([mol.GetNumConformers() >= 1 for mol in mols]):
        raise ValueError("One or more input molecules is missing a conformer.")

    # Make a copy of the molecules since they are going to be modified
    if copy:
        mols = [dm_mol.copy_mol(mol) for mol in mols]

    # Split ref and probe mols
    mol_ref = mols[ref_id]
    mol_probes = mols

    if backend == "crippenO3A":

        # Compute Crippen contributions for every atoms and molecules
        crippen_contribs = [rdMolDescriptors._CalcCrippenContribs(mol) for mol in mol_probes]
        crippen_contrib_ref = crippen_contribs[ref_id]
        crippen_contrib_probes = crippen_contribs

        # Loop and align
        # NOTE(hadim): we could eventually parallelize this if that's needed.

        scores = []
        for i, mol in enumerate(mol_probes):

            crippenO3A = rdMolAlign.GetCrippenO3A(
                prbMol=mol,
                refMol=mol_ref,
                prbCrippenContribs=crippen_contrib_probes[i],
                refCrippenContribs=crippen_contrib_ref,
                prbCid=conformer_id,
                refCid=conformer_id,
                maxIters=50,
            )
            crippenO3A.Align()

            scores.append(crippenO3A.Score())

    elif backend == "O3A":

        # Add hydrogens first
        mol_probes = [dm_mol.add_hs(mol, add_coords=True) for mol in mol_probes]
        mol_ref = dm_mol.add_hs(mol_ref, add_coords=True)

        # Compute MMFF params for every molecules
        mmff_params = [rdForceFieldHelpers.MMFFGetMoleculeProperties(mol) for mol in mol_probes]

        # Split reference and probe molecules
        mmff_params_ref = mmff_params[ref_id]
        mmff_params_probes = mmff_params

        # Loop and align
        # NOTE(hadim): we could eventually parallelize this if that's needed.

        scores = []
        for i, mol in enumerate(mol_probes):

            pyO3A = rdMolAlign.GetO3A(
                prbMol=mol,
                refMol=mol_ref,
                prbPyMMFFMolProperties=mmff_params_probes[i],
                refPyMMFFMolProperties=mmff_params_ref,
                prbCid=conformer_id,
                refCid=conformer_id,
                maxIters=50,
            )
            pyO3A.Align()

            scores.append(pyO3A.Score())

        # Remove the hydrogens
        mol_probes = [dm_mol.remove_hs(mol) for mol in mol_probes]

    else:
        raise ValueError(f"Backend {backend} not supported.")

    scores = np.array(scores)

    return mol_probes, scores
