import sys
import os

from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import RDConfig
from rdkit.Chem import Lipinski
from rdkit.Chem import rdmolops
from rdkit.Chem import Crippen


from .. import Mol
from ..convert import from_smarts


def _sasscorer(mol: Mol):

    sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
    try:
        import sascorer  # type:ignore
    except ImportError:
        raise ImportError(
            "Could not import sascorer. If you installed rdkit-pypi with `pip`, please uninstall it and reinstall rdkit with `conda` or `mamba`."
        )

    return sascorer.calculateScore(mol)


_AROMATIC_QUERY = from_smarts("a")

mw = rdMolDescriptors.CalcExactMolWt
fsp3 = rdMolDescriptors.CalcFractionCSP3
tpsa = rdMolDescriptors.CalcTPSA
qed = Descriptors.qed
clogp = Descriptors.MolLogP  # type: ignore
sas = _sasscorer
formal_charge = rdmolops.GetFormalCharge
refractivity = Crippen.MolMR

n_hba = rdMolDescriptors.CalcNumHBA
n_hbd = rdMolDescriptors.CalcNumHBD
n_lipinski_hba = rdMolDescriptors.CalcNumLipinskiHBA
n_lipinski_hbd = rdMolDescriptors.CalcNumLipinskiHBD
n_rings = rdMolDescriptors.CalcNumRings
n_hetero_atoms = rdMolDescriptors.CalcNumHeteroatoms
n_heavy_atoms = Descriptors.HeavyAtomCount  # type: ignore
n_rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds
n_radical_electrons = Descriptors.NumRadicalElectrons
n_NHOH = Lipinski.NHOHCount
n_NO = Lipinski.NOCount

n_aliphatic_carbocycles = Lipinski.NumAliphaticCarbocycles  # type: ignore
n_aliphatic_heterocyles = Lipinski.NumAliphaticHeterocycles  # type: ignore
n_aliphatic_rings = Lipinski.NumAliphaticRings  # type: ignore

n_aromatic_carbocycles = Lipinski.NumAromaticCarbocycles  # type: ignore
n_aromatic_heterocyles = Lipinski.NumAromaticHeterocycles  # type: ignore
n_aromatic_rings = Lipinski.NumAromaticRings  # type: ignore

n_saturated_carbocycles = Lipinski.NumSaturatedCarbocycles  # type: ignore
n_saturated_heterocyles = Lipinski.NumSaturatedHeterocycles  # type: ignore
n_saturated_rings = Lipinski.NumSaturatedRings  # type: ignore


def n_rigid_bonds(mol: Mol) -> int:
    """Compute the number of rigid bonds in a molecule.

    Rigid bonds are bonds that are not single and not in rings.

    Args:
        mol: A molecule.

    Returns:
        n_rigid_bonds: number of rigid bonds in the molecule
    """
    non_rigid_bonds_count = from_smarts("*-&!@*")
    n_rigid_bonds = mol.GetNumBonds() - len(mol.GetSubstructMatches(non_rigid_bonds_count))
    return n_rigid_bonds


def n_aromatic_atoms(mol: Mol) -> int:
    """Calculate the number of aromatic atoms."""
    matches = mol.GetSubstructMatches(_AROMATIC_QUERY)
    return len(matches)


def n_aromatic_atoms_proportion(mol: Mol) -> int:
    """Calculate the aromatic proportion: # aromatic atoms/#atoms total.

    Args:
        mol: A molecule.

    Only heavy atoms are considered.
    """
    return n_aromatic_atoms(mol) / mol.GetNumHeavyAtoms()


def n_stereo_centers(mol: Mol) -> int:
    """Compute the number of stereocenters in a molecule.

    Args:
        mol: A molecule.

    Returns:
        n_stero_center: number of stereocenters in the molecule
    """
    n = 0
    try:
        rdmolops.FindPotentialStereo(mol, cleanIt=False)
        n = rdMolDescriptors.CalcNumAtomStereoCenters(mol)
    except:
        pass
    return n


def n_charged_atoms(mol: Mol) -> int:
    """Compute the number of charged atoms in a molecule.

    Args:
        mol: A molecule.

    Returns:
        n_charged_atoms: number of charged atoms in the molecule
    """
    return sum([at.GetFormalCharge() != 0 for at in mol.GetAtoms()])
