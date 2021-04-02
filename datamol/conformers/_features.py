from typing import Union
from typing import List

import numpy as np

from rdkit import Chem

import datamol as dm


@dm.utils.decorators.disable_on_os("win")
def sasa(
    mol: Chem.rdchem.Mol,
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
        mol: the molecule with the conformers.
    """
    from rdkit.Chem import rdFreeSASA

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


def get_coords(mol: Chem.rdchem.Mol, conf_id: int = -1):
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
    mol: Chem.rdchem.Mol,
    use_atoms: bool = True,
    digits: int = None,
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
    coords = get_coords(mol)
    atom_weight = np.ones((coords.shape[0]))

    if use_atoms:
        atom_weight = np.array([atom.GetMass() for atom in mol.GetAtoms()])

    atom_weight = atom_weight[:, None]
    atom_weight /= atom_weight.sum()
    center = (coords * atom_weight).sum(axis=0)

    if digits is not None:
        center = center.round(digits)

    return center
