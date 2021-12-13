from rdkit import Chem
from rdkit.Chem import rdmolops

from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
from rdkit.Chem.EnumerateStereoisomers import StereoEnumerationOptions

import datamol as dm

from ._structural import IsomerEnumerator


def enumerate_tautomers(mol: Chem.rdchem.Mol, n_variants: int = 20):
    """Enumerate the possible tautomers of the current molecule.

    Args:
        mol: The molecule whose state we should enumerate.
        n_variants: The maximum amount of molecules that should be returned.
    """
    enumerator = rdMolStandardize.TautomerEnumerator()
    enumerator.SetMaxTautomers(n_variants)
    tautomers = enumerator.Enumerate(mol)
    return list(tautomers)


def enumerate_stereoisomers(
    mol,
    n_variants: int = 20,
    undefined_only: bool = False,
    rationalise: bool = True,
):
    """Enumerate the stereocenters and bonds of the current molecule.

    Original source: the `openff-toolkit` lib.

    Warning: this function can be computationnaly intensive.

    Args:
        mol: The molecule whose state we should enumerate.
        n_variants: The maximum amount of molecules that should be returned.
        undefined_only: If we should enumerate all stereocenters and bonds or only those
            with undefined stereochemistry.
        rationalise: If we should try to build and rationalise the molecule to ensure it
            can exist.
    """

    # safety first
    mol = dm.copy_mol(mol)

    # in case any bonds/centers are missing stereo chem flag it here
    Chem.AssignStereochemistry(mol, force=False, flagPossibleStereoCenters=True, cleanIt=True)  # type: ignore
    Chem.FindPotentialStereoBonds(mol, cleanIt=True)  # type: ignore

    # set up the options
    stereo_opts = StereoEnumerationOptions(
        tryEmbedding=rationalise,
        onlyUnassigned=undefined_only,
        maxIsomers=n_variants,
    )

    try:
        isomers = tuple(EnumerateStereoisomers(mol, options=stereo_opts))
    except:
        # NOTE(hadim): often got "Stereo atoms should be specified before specifying CIS/TRANS bond stereochemistry"
        # for the ligand of reference (coming from the PDB). Not sure how to handle that.
        isomers = []

    variants = []
    for isomer in isomers:
        # isomer has CIS/TRANS tags so convert back to E/Z
        Chem.SetDoubleBondNeighborDirections(isomer)  # type: ignore
        Chem.AssignStereochemistry(isomer, force=True, cleanIt=True)  # type: ignore
        variants.append(isomer)

    return variants


def enumerate_structisomers(
    mol,
    n_variants: int = 20,
    allow_cycle: bool = False,
    allow_double_bond: bool = False,
    allow_triple_bond: bool = False,
    depth: int = None,
):
    """Enumerate the structural isomers of the input molecule

    Warning: this function can be computationnaly intensive.

    Args:
        mol: The molecule whose state we should enumerate.
        n_variants: The maximum amount of molecules that should be returned.
        allow_cycle: whether to allow transformation involving cycle creation
        allow_double_bond: whether to allow transformation involving at least one double bond
        allow_triple_bond: whether to allow transformation involving at least one triple bond
        depth: depth of the search, use a sensible value to decrease search time. Will fall back to number of atoms
    """
    mol = dm.copy_mol(mol)

    enumerator = IsomerEnumerator(
        allow_cycle=allow_cycle,
        allow_double_bond=allow_double_bond,
        allow_triple_bond=allow_triple_bond,
    )
    if depth is None:
        depth = mol.GetNumAtoms()
    variants = enumerator(mol, max_mols=n_variants, include_input=False, depth=depth)
    return [dm.to_mol(x) for x in variants]


def remove_stereochemistry(mol: dm.Mol, copy: bool = True):
    """Removes all stereochemistry info from the molecule.

    Args:
        mol: a molecule
        copy: Whether to copy the molecule.
    """

    if copy is True:
        mol = dm.copy_mol(mol)
    rdmolops.RemoveStereochemistry(mol)
    return mol


def canonical_tautomer(mol: dm.Mol):
    """Get the canonical tautomer of the current molecule.

    Args:
        mol: A molecule.
    """
    enumerator = rdMolStandardize.TautomerEnumerator()
    return enumerator.Canonicalize(mol)
