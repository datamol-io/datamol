from typing import List

from rdkit.Chem import rdFMCS

import datamol as dm


def find_mcs_with_details(
    mols: List[dm.Mol],
    maximize_bonds: bool = True,
    threshold: float = 1.0,
    timeout: int = 5,
    verbose: bool = False,
    match_valences: bool = False,
    ring_matches_ring_only: bool = True,
    complete_rings_only: bool = False,
    match_chiral_tag: bool = False,
    seed_smarts: str = "",
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
        kwargs: Additional arguments for the MCS.
    """

    # NOTE(hadim): eventually propagate enums from
    # AtomCompare and BondCompare. Details at
    # https://www.rdkit.org/docs/source/rdkit.Chem.rdFMCS.html

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
