from typing import Union
from typing import List
from typing import Dict
from typing import Any

import copy

from loguru import logger

import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdFreeSASA
from rdkit.Chem import TorsionFingerprints
from rdkit.Chem import Descriptors
from rdkit.ML.Cluster import Butina

import datamol as dm


def generate(
    mol: Chem.Mol,
    n_confs: int = None,
    method: str = None,
    align_conformers: bool = True,
    minimize_energy: bool = True,
    energy_iterations: int = 500,
    warning_not_converged: int = 10,
    random_seed: int = 19,
    verbose: bool = False,
) -> Chem.Mol:
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
            number of rotatable bonds by default. Defaults to None.
        method: RDKit method to use for embedding. Choose among
            ["ETDG", "ETKDG", "ETKDGv2", "ETKDGv3"]. If None, "ETKDGv3" is used. Default to None.
        align_conformers: Wehther to align conformer. Note that this is done
            BEFORE the energy minimization procedure. Defaults to True.
        minimize_energy: Wether to minimize conformer's energies using UFF.
            Disable to generate conformers much faster. Defaults to True.
        energy_iterations: Maximum number of iterations during the energy minimization procedure.
            It corresponds to the `maxIters` argument in RDKit. Defaults to 500.
        warning_not_converged: Wether to log a warning when the number of not converged conformers
            during the minimization is higher than `warning_not_converged`. Only works when `verbose` is set to True. Disable with 0. Defaults to 10.
        random_seed: Set to None or -1 to disable. Defaults to 19.
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

    # It's probably best to compute conformers using the hydrogen
    # atoms. Even if rdkit does not always use those.
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

    # Align conformers to each others
    if align_conformers:
        Chem.rdMolAlign.AlignMolConformers(mol)

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
        ordered_conformers = [conf for _, conf in sorted(zip(energies, mol.GetConformers()))]
        mol.RemoveAllConformers()
        [mol.AddConformer(conf, assignId=True) for conf in ordered_conformers]

    return mol


def cluster(
    mol: Chem.Mol,
    method: str = None,
    distance_threshold: float = 1,
    return_centroids: bool = True,
    method_kwargs: Dict[Any, Any] = None,
) -> Union[List[Chem.Mol], Chem.Mol]:
    """Cluster molecule's conformers.

    Args:
        mol: a molecule
        method: Distance method to use from ["RMS", "TFD"].
            Use None for "RMS", Defaults to None.
            - TFD uses `TorsionFingerprints.GetTFDMatrix`.
            - RMS uses `AllChem.GetConformerRMSMatrix`.
        distance_threshold: Threshold for the clustering. Defaults to 1.
        return_centroids: If True, return one molecule with centroid conformers
            only. If False return a list of molecules per cluster with all the conformers of the cluster. Defaults to True.
        method_kwargs: Additional method to pass to the clustering function.
    """

    AVAILABLE_METHODS = ["RMS", "TFD"]

    # Clone molecule
    mol = copy.deepcopy(mol)

    if method_kwargs is None:
        method_kwargs = {}

    if method is None:
        method = "RMS"

    if method == "TFD":
        kwargs = {}
        kwargs["useWeights"] = False
        kwargs["maxDev"] = "equal"
        kwargs["symmRadius"] = 2
        kwargs["ignoreColinearBonds"] = True
        kwargs.update(method_kwargs)
        dmat = TorsionFingerprints.GetTFDMatrix(mol, **kwargs)
    elif method == "RMS":
        kwargs = {}
        kwargs["prealigned"] = False
        kwargs.update(method_kwargs)
        dmat = AllChem.GetConformerRMSMatrix(mol, prealigned=False)
    else:
        raise ValueError(f"The method {method} is not supported. Use from {AVAILABLE_METHODS}")

    conf_clusters = Butina.ClusterData(
        dmat,
        nPts=mol.GetNumConformers(),
        distThresh=distance_threshold,
        isDistData=True,
        reordering=False,
    )

    # Collect centroid of each cluster (first element of the list)
    centroids = [indices[0] for indices in conf_clusters]

    if return_centroids:
        # Keep only centroid conformers
        confs = [mol.GetConformers()[i] for i in centroids]
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


def sasa(
    mol: Chem.Mol,
    conf_id: Union[int, List[int]] = None,
    n_jobs: int = 1,
) -> np.ndarray:
    """Compute Solvent Accessible Surface Area of all the conformers
    using FreeSASA (https://freesasa.github.io/). Values are returned
    as an array and also stored within each conformer as a property
    called `rdkit_free_sasa`.

    Example:

    ```python
    smiles = "O=C(C)Oc1ccccc1C(=O)O"
    mol = dm.to_mol(smiles)
    mol = dm.conformers.generate(mol)

    # Compute SASA for all the conformers without parallelization
    sasa_values = dm.conformers.sasa(mol, conf_id=None, n_jobs=1)

    # If minimization has been enabled (default to True)
    # you can access the computed energy.
    conf = mol.GetConformer(0)
    props = conf.GetPropsAsDict()
    print(props)
    # {'rdkit_uff_energy': 1.7649408317784008}
    ```

    Args:
        mol: a molecule
        conf_id: Id of the conformers to compute. If None, compute all.
        n_jobs: Number of jobs for parallelization. Set to 1 to disable
            and -1 to use all cores.

    Returns:
        mol (Chem.Mol): the molecule with the conformers.
    """

    if mol.GetNumConformers() == 0:
        raise ValueError(
            "The molecule has 0 conformers. You can generate conformers with `dm.conformers.generate(mol)`."
        )

    # Get Van der Waals radii (angstrom)
    radii = [dm.PERIODIC_TABLE.GetRvdw(atom.GetAtomicNum()) for atom in mol.GetAtoms()]

    # Which conformers to compute
    conf_ids = []
    if conf_id is None:
        # If None compute for all the conformers
        conf_ids = list(range(mol.GetNumConformers()))  # type: ignore
    elif isinstance(conf_id, int):
        conf_ids = [conf_id]
    else:
        conf_ids = conf_id

    # Compute solvent accessible surface area
    def _get_sasa(i):
        conf = mol.GetConformer(i)
        sasa = rdFreeSASA.CalcSASA(mol, radii, confIdx=conf.GetId())
        conf.SetDoubleProp("rdkit_free_sasa", sasa)
        return sasa

    runner = dm.JobRunner(n_jobs=n_jobs)
    sasa_values = runner(_get_sasa, conf_ids)
    return np.array(sasa_values)


def rmsd(mol: Chem.Mol) -> np.ndarray:
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
            rmsd = Chem.rdMolAlign.AlignMol(prbMol=mol, refMol=mol, prbCid=i, refCid=j)
            rmsds.append(rmsd)
    return np.array(rmsds).reshape(n_confs, n_confs)
