import numpy as np

from rdkit import Chem


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
