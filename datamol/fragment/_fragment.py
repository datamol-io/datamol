from typing import Set
from typing import Optional
from typing import Any

from rdkit import Chem
from rdkit.Chem import BRICS
from rdkit.Chem import Recap
from rdkit.Chem import rdMMPA

from rdkit.Chem.Fraggle import FraggleSim

import datamol as dm


def brics(
    mol: Chem.rdchem.Mol,
    singlepass: bool = True,
    remove_parent: bool = False,
    sanitize: bool = True,
    fix: bool = True,
):
    """Run BRICS on the molecules and potentially fix dummy atoms.

    Args:
        mol: a molecule.
        singlepass: Single pass for `BRICSDecompose`.
        remove_parent: Remove parent from the fragments.
        sanitize: Wether to sanitize the fragments.
        fix: Wether to fix the fragments.
    """
    frags = BRICS.BRICSDecompose(mol, returnMols=True, singlePass=singlepass)
    frags = list(frags)

    if fix:
        frags = [dm.fix_mol(x) for x in frags]
    if sanitize:
        frags = [dm.sanitize_mol(x) for x in frags]
    if remove_parent:
        frags.pop(0)

    frags = [x for x in frags if x is not None]

    return frags


def frag(
    mol: Chem.rdchem.Mol,
    remove_parent: bool = False,
    sanitize: bool = True,
    fix: bool = True,
):
    """Generate all possible fragmentation of a molecule.

    Args:
        mol: a molecule.
        remove_parent: Remove parent from the fragments.
        sanitize: Wether to sanitize the fragments.
        fix: Wether to fix the fragments.
    """
    frags = FraggleSim.generate_fraggle_fragmentation(mol)

    smiles = set([])
    for seq in frags:
        smiles |= {s.strip() for s in seq.split(".")}

    smiles = list(sorted(smiles, reverse=True))
    frags = [dm.to_mol(s) for s in smiles]

    if fix:
        frags = [dm.fix_mol(x) for x in frags]
    if sanitize:
        frags = [dm.sanitize_mol(x) for x in frags]

    frags = [x for x in frags if x is not None]

    if remove_parent:
        return frags
    return [mol] + frags


def recap(
    mol: Chem.rdchem.Mol,
    remove_parent: bool = False,
    sanitize: bool = True,
    fix: bool = True,
):
    """Fragment the molecule using the recap algorithm.

    Args:
        mol: a molecule.
        remove_parent: Remove parent from the fragments.
        sanitize: Wether to sanitize the fragments.
        fix: Wether to fix the fragments.
    """
    res = Recap.RecapDecompose(mol)
    frags = [dm.to_mol(x) for x in res.GetAllChildren().keys()]

    if fix:
        frags = [dm.fix_mol(x) for x in frags]
    if sanitize:
        frags = [dm.sanitize_mol(x) for x in frags]

    frags = [x for x in frags if x is not None]

    if remove_parent:
        return frags
    return [mol] + frags


def anybreak(
    mol: Chem.rdchem.Mol,
    remove_parent: bool = False,
    sanitize: bool = True,
    fix: bool = True,
):
    """Fragment molecule by applying brics first, then fall back to frag.

    Args:
        mol: a molecule.
        remove_parent: Remove parent from the fragments.
        sanitize: Wether to sanitize the fragments.
        fix: Wether to fix the fragments.
    """
    frags = []
    try:
        frags = brics(mol, fix=fix, remove_parent=remove_parent, sanitize=sanitize)
    except:
        pass

    if len(frags) == 0:
        frags = frag(mol, remove_parent=remove_parent, sanitize=sanitize, fix=fix)

    return frags


def mmpa_frag(
    mol: dm.Mol,
    pattern: Optional[str] = None,
    max_cut: int = 1,
    max_bond_cut: int = 20,
    h_split: bool = False,
) -> Optional[Set[dm.Mol]]:
    """Fragment molecule on specific bonds suitable for a MMPA analysis.

    Args:
        mol: Molecule to fragment.
        pattern: Bond pattern to split on. Will use default rdkit pattern
            '[#6+0;!$(*=,#[!#6])]!@!=!#[*]' if not provided.
        max_cut: Number of cuts.
        max_bond_cut: Maximum number of bond to cut. Default to 20.
        h_split:  Whether to split at hydrogen position too.
            This is equivalent to enabling the addition of new fragments.

    Returns:
        List of fragments.
    """

    frags = []
    if pattern is None:
        frags = rdMMPA.FragmentMol(
            mol,
            maxCuts=max_cut,
            resultsAsMols=False,
            maxCutBonds=max_bond_cut,
        )
    elif pattern:
        frags = rdMMPA.FragmentMol(
            mol,
            pattern=pattern,
            maxCuts=max_cut,
            resultsAsMols=False,
            maxCutBonds=max_bond_cut,
        )

    if h_split:
        mol = dm.add_hs(mol)
        frags += rdMMPA.FragmentMol(
            mol,
            pattern="[#1]!@!=!#[!#1]",
            maxCuts=1,
            resultsAsMols=False,
            maxCutBonds=max_bond_cut,
        )
    return set(frags)


def mmpa_cut(mol: dm.Mol, rdkit_pattern: bool = False) -> Optional[Set[Any]]:
    """Cut molecules to perform mmpa analysis later

    Args:
        mol: Molecule to fragment.
        rdkit_pattern: Whether to perform the fragmentation
            using the default rdkit pattern: [#6+0;!$(*=, #[!#6])]!@!=!#[*]"

    Returns:
        List of 'smiles,core,chains'
    """

    if mol is None:
        return mol

    outlines = set()

    smiles = dm.to_smiles(mol)

    if rdkit_pattern:
        frags = mmpa_frag(mol, max_cut=3, max_bond_cut=30)
    else:
        # heavy atoms
        frags = mmpa_frag(mol, pattern="[!#1]!@!=!#[!#1]", max_cut=4, max_bond_cut=30)
        frags.update(mmpa_frag(mol, pattern="[!#1]!@!=!#[!#1]", max_cut=3, max_bond_cut=30))

    frags = set(frags)
    for core, chains in frags:
        output = f"{smiles},{core},{chains}\n"
        outlines.add(output)

    # hydrogen splitting
    mol = dm.add_hs(mol)
    smiles = dm.to_smiles(mol)

    n = mol.GetNumHeavyAtoms()
    if n < 60:
        frags = mmpa_frag(mol, pattern=None, max_cut=1, max_bond_cut=100, h_split=True)
        for core, chains in frags:
            output = f"{smiles},{core},{chains}\n"
            outlines.add(output)

    return outlines
