from typing import List

from rdkit.Chem import rdFMCS

import datamol as dm

ALLOWED_ATOM_COMPARE = ["CompareAny", "CompareAnyHeavyAtom", "CompareElements", "CompareIsotopes"]
ALLOWED_BOND_COMPARE = ["CompareAny", "CompareOrder", "CompareOrderExact"]
ALLOWED_RING_COMPARE = ["IgnoreRingFusion", "PermissiveRingFusion", "StrictRingFusion"]


def find_mcs_with_details(
    mols: List[dm.Mol],
    maximize_bonds: bool = True,
    threshold: float = 0.0,
    timeout: int = 5,
    verbose: bool = False,
    match_valences: bool = False,
    ring_matches_ring_only: bool = True,
    complete_rings_only: bool = False,
    match_chiral_tag: bool = False,
    seed_smarts: str = "",
    atom_compare: str = "CompareElements",
    bond_compare: str = "CompareOrder",
    ring_compare: str = "IgnoreRingFusion",
    **kwargs,
):
    """Find the maximum common substructure from a list of molecules.

    Args:
        mols: List of molecules.
        maximize_bonds: Maximize the number of bonds in the substructure.
        threshold: The threshold for the MCS (between 0 and 1).
        timeout: The timeout for the MCS.
        verbose: Whether to enable verbose mode.
        match_valences: Whether to match valences.
        ring_matches_ring_only: Whether to match rings only.
        complete_rings_only: Whether to match complete rings only.
        match_chiral_tag: Whether to match chiral tags.
        seed_smarts: The seed SMARTS.
        atom_compare: One of "CompareAny", "CompareAnyHeavyAtom", "CompareElements",
            "CompareIsotopes".
        bond_compare: One of "CompareAny", "CompareOrder", "CompareOrderExact".
        ring_compare: One of "IgnoreRingFusion", "PermissiveRingFusion", "StrictRingFusion".
        kwargs: Additional arguments for the MCS.
    """

    if atom_compare not in ALLOWED_ATOM_COMPARE:
        raise ValueError(f"atom_compare must be one of {ALLOWED_ATOM_COMPARE}")

    if bond_compare not in ALLOWED_BOND_COMPARE:
        raise ValueError(f"bond_compare must be one of {ALLOWED_BOND_COMPARE}")

    if ring_compare not in ALLOWED_RING_COMPARE:
        raise ValueError(f"ring_compare must be one of {ALLOWED_RING_COMPARE}")

    args = {}
    args["maximizeBonds"] = maximize_bonds
    args["threshold"] = threshold
    args["timeout"] = timeout
    args["verbose"] = verbose
    args["matchValences"] = match_valences
    args["ringMatchesRingOnly"] = ring_matches_ring_only
    args["completeRingsOnly"] = complete_rings_only
    args["matchChiralTag"] = match_chiral_tag
    args["seedSmarts"] = seed_smarts
    args["atomCompare"] = rdFMCS.AtomCompare.names[atom_compare]
    args["bondCompare"] = rdFMCS.BondCompare.names[bond_compare]
    args["ringCompare"] = rdFMCS.RingCompare.names[ring_compare]

    args.update(kwargs)

    mcs = rdFMCS.FindMCS(mols, **args)
    return mcs


def find_mcs(
    mols: List[dm.Mol],
    maximize_bonds: bool = True,
    threshold: float = 1.0,
    timeout: int = 5,
    verbose: bool = False,
    match_valences: bool = False,
    ring_matches_ring_only: bool = False,
    complete_rings_only: bool = False,
    match_chiral_tag: bool = False,
    seed_smarts: str = "",
    **kwargs,
):
    """Find the maximum common substructure from a list of molecules.

    Args:
        mols: List of molecules.
        maximize_bonds: Maximize the number of bonds in the substructure.
        threshold: The threshold for the MCS.
        timeout: The timeout for the MCS.
        verbose: Whether to enable verbose mode.
        match_valences: Whether to match valences.
        ring_matches_ring_only: Whether to match rings only.
        complete_rings_only: Whether to match complete rings only.
        match_chiral_tag: Whether to match chiral tags.
        seed_smarts: The seed SMARTS.
        kwargs: Additional arguments for the MCS.
    """

    mcs = find_mcs_with_details(
        mols=mols,
        maximize_bonds=maximize_bonds,
        threshold=threshold,
        timeout=timeout,
        verbose=verbose,
        match_valences=match_valences,
        ring_matches_ring_only=ring_matches_ring_only,
        complete_rings_only=complete_rings_only,
        match_chiral_tag=match_chiral_tag,
        seed_smarts=seed_smarts,
        **kwargs,
    )
    return mcs.smartsString
