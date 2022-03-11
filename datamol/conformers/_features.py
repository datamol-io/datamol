from typing import Optional, Union
from typing import List
from typing import Optional

import numpy as np

from rdkit import Chem

from ..types import Mol
from ..utils.jobs import JobRunner
from ..utils import decorators
from ..mol import PERIODIC_TABLE
from ..mol import copy_mol


@decorators.disable_on_os("win")
def sasa(
    mol: Mol,
    conf_id: Optional[Union[int, List[int]]] = None,
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
        mol: the molecule with the conformers.
    """
    from rdkit.Chem import rdFreeSASA

    if mol.GetNumConformers() == 0:
        raise ValueError(
            "The molecule has 0 conformers. You can generate conformers with `dm.conformers.generate(mol)`."
        )

    # Get Van der Waals radii (angstrom)
    radii = [PERIODIC_TABLE.GetRvdw(atom.GetAtomicNum()) for atom in mol.GetAtoms()]

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

    runner = JobRunner(n_jobs=n_jobs)
    sasa_values = runner(_get_sasa, conf_ids)
    return np.array(sasa_values)


def get_coords(mol: Mol, conf_id: int = -1):
    """Get the coordinate of a conformer of a molecule.

    Args:
        mol: a molecule.
        conf_id: a conformer id.
    """

    if mol.GetNumConformers() == 0:
        raise ValueError("Molecule does not have any conformers.")

    conf = mol.GetConformer(id=conf_id)
    return conf.GetPositions()


def center_of_mass(
    mol: Mol,
    use_atoms: bool = True,
    digits: Optional[int] = None,
    conf_id: int = -1,
) -> np.ndarray:
    """Compute the center of mass of a conformer of a molecule.

    Args:
        mol: a molecule
        use_atoms: Whether to compute the true center of mass or the geometrical center.
        digits: Number of digits to round to.
        conf_id: the conformer id.

    Returns
        cm: Center of mass or geometrical center
    """
    coords = get_coords(mol, conf_id=conf_id)
    atom_weight = np.ones((coords.shape[0]))

    if use_atoms:
        atom_weight = np.array([atom.GetMass() for atom in mol.GetAtoms()])

    atom_weight = atom_weight[:, None]
    atom_weight /= atom_weight.sum()
    center = (coords * atom_weight).sum(axis=0)

    if digits is not None:
        center = center.round(digits)

    return center


def keep_conformers(
    mol: Mol,
    indices_to_keep: Union[int, List[int]] = -1,
    assign_id: bool = True,
    copy: bool = True,
):
    """Keep on the specified conformer(s) in `indices_to_keep`.

    Args:
        mol: A molecule.
        indices_to_keep: A indice or a least of indices of conformers to keep.
        assign_id: Whether to assign the kept conformers an id or keep the original one.
        copy: Whether to copy the molecule or not.
    """

    if copy:
        mol = copy_mol(mol)

    if not isinstance(indices_to_keep, list):
        indices_to_keep = [indices_to_keep]

    # Extract conformers to keep
    confs_to_keep = [mol.GetConformer(conf_id) for conf_id in indices_to_keep]

    # Copy current mol and remove all conformers
    mol2 = copy_mol(mol)
    mol2.RemoveAllConformers()

    # Add conformers
    _ = [mol2.AddConformer(conf, assignId=assign_id) for conf in confs_to_keep]

    # Cleanup
    mol = mol2

    return mol
