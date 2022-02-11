import time

from rdkit import Chem
from rdkit.Chem import rdmolops

from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
from rdkit.Chem.EnumerateStereoisomers import StereoEnumerationOptions

import datamol as dm

from ._structural import IsomerEnumerator


def enumerate_tautomers(
    mol: dm.Mol,
    n_variants: int = 20,
    max_transforms: int = 1000,
    reassign_stereo: bool = True,
    remove_bond_stereo: bool = True,
    remove_sp3_stereo: bool = True,
):
    """Enumerate the possible tautomers of the current molecule.

    Args:
        mol: The molecule whose state we should enumerate.
        n_variants: The maximum amount of molecules that should be returned.
        max_transforms: Set the maximum number of transformations to be applied. This limit is usually
            hit earlier than the n_variants limit and leads to a more linear scaling
            of CPU time with increasing number of tautomeric centers (see Sitzmann et al.).
        reassign_stereo: Whether to reassign stereo centers.
        remove_bond_stereo: Whether stereochemistry information is removed from double bonds involved in tautomerism.
            This means that enols will lose their E/Z stereochemistry after going through tautomer enumeration because
            of the keto-enolic tautomerism.
        remove_sp3_stereo: Whether stereochemistry information is removed from sp3 atoms involved in tautomerism. This
            means that S-aminoacids will lose their stereochemistry after going through tautomer enumeration because
            of the amido-imidol tautomerism.
    """
    enumerator = rdMolStandardize.TautomerEnumerator()

    # Configure the enumeration
    enumerator.SetMaxTautomers(n_variants)
    enumerator.SetMaxTransforms(max_transforms)
    enumerator.SetReassignStereo(reassign_stereo)
    enumerator.SetRemoveBondStereo(remove_bond_stereo)
    enumerator.SetRemoveSp3Stereo(remove_sp3_stereo)

    tautomers = enumerator.Enumerate(mol)
    return list(tautomers)


def enumerate_stereoisomers(
    mol: dm.Mol,
    n_variants: int = 20,
    undefined_only: bool = False,
    rationalise: bool = True,
    timeout_seconds: int = None,
    clean_it: bool = True,
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
        timeout_seconds: The maximum amount of time to spend on enumeration. None
            will disable the timeout. Note that the timeout might be inaccurate as a running single variant
            computation is not stopped when the duration is reached.
        clean_it: A flag for assigning stereochemistry. If True, it will remove previous stereochemistry
            markings on the bonds.
    """

    # safety first
    mol = dm.copy_mol(mol)

    # in case any bonds/centers are missing stereo chem flag it here
    Chem.AssignStereochemistry(mol, force=False, flagPossibleStereoCenters=True, cleanIt=clean_it)  # type: ignore
    Chem.FindPotentialStereoBonds(mol, cleanIt=clean_it)  # type: ignore

    # set up the options
    stereo_opts = StereoEnumerationOptions(
        tryEmbedding=rationalise,
        onlyUnassigned=undefined_only,
        maxIsomers=n_variants,
        unique=True,
    )

    isomers_iterator = EnumerateStereoisomers(mol, options=stereo_opts)

    start = time.time()
    duration = 0

    isomers = []
    while timeout_seconds is None or duration < timeout_seconds:

        try:
            isomer = next(isomers_iterator)
            isomers.append(isomer)
        except StopIteration:
            break

        duration = time.time() - start

    variants = []
    for isomer in isomers:
        # isomer has CIS/TRANS tags so convert back to E/Z
        Chem.SetDoubleBondNeighborDirections(isomer)  # type: ignore
        Chem.AssignStereochemistry(isomer, force=True, cleanIt=clean_it)  # type: ignore
        variants.append(isomer)

    return variants


def enumerate_structisomers(
    mol: dm.Mol,
    n_variants: int = 20,
    allow_cycle: bool = False,
    allow_double_bond: bool = False,
    allow_triple_bond: bool = False,
    depth: int = None,
    timeout_seconds: int = None,
):
    """Enumerate the structural isomers of the input molecule

    Warning: this function can be computationnaly intensive.

    Args:
        mol: The molecule whose state we should enumerate.
        n_variants: The maximum amount of molecules that should be returned.
        allow_cycle: whether to allow transformation involving cycle creation.
        allow_double_bond: whether to allow transformation involving at least one double bond.
        allow_triple_bond: whether to allow transformation involving at least one triple bond.
        depth: depth of the search, use a sensible value to decrease search time. Will fall back to number of atoms.
        timeout_seconds: The maximum amount of time to spend on enumeration. None
            will disable the timeout. Note that the timeout might be inaccurate as a running single variant
            computation is not stopped when the duration is reached.
    """

    mol = dm.copy_mol(mol)

    enumerator = IsomerEnumerator(
        allow_cycle=allow_cycle,
        allow_double_bond=allow_double_bond,
        allow_triple_bond=allow_triple_bond,
    )

    if depth is None:
        depth = mol.GetNumAtoms()

    isomers_iterator = enumerator(
        mol,
        max_mols=n_variants,
        include_input=False,
        depth=depth,
    )

    start = time.time()
    duration = 0

    isomers = []
    while timeout_seconds is None or duration < timeout_seconds:

        try:
            isomer = next(isomers_iterator)
            isomers.append(dm.to_mol(isomer))
        except StopIteration:
            break

        duration = time.time() - start

    return isomers


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
