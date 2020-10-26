from rdkit import Chem
from rdkit.Chem import BRICS
from rdkit.Chem import Recap
from rdkit.Chem import rdRGroupDecomposition
from rdkit.Chem import rdMMPA

from rdkit.Chem.Fraggle import FraggleSim

import datamol as dm


def brics(
    mol: Chem.Mol,
    singlepass: bool = True,
    remove_parent: bool = False,
    sanitize: bool = True,
    fix: bool = True,
):
    """Run BRICS on the molecules and potentially fix dummy atoms.

    Args:
        mols (Chem.Mol): a molecule.
        singlepass (bool, optional): Single pass for `BRICSDecompose`.
        remove_parent (bool, optional): Remove parent from the fragments.
        sanitize (bool, optional): Wether to sanitize the fragments.
        fix (bool, optional): Wether to fix the fragments.
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
    mol: Chem.Mol,
    remove_parent: bool = False,
    sanitize: bool = True,
    fix: bool = True,
):
    """Fragment molecule using FraggleSim on all bonds.

    Args:
        mols (Chem.Mol): a molecule.
        remove_parent (bool, optional): Remove parent from the fragments.
        sanitize (bool, optional): Wether to sanitize the fragments.
        fix (bool, optional): Wether to fix the fragments.
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
    mol: Chem.Mol,
    remove_parent: bool = False,
    sanitize: bool = True,
    fix: bool = True,
):
    """Fragment the molecule using the recap algorithm.

    Args:
        mols (Chem.Mol): a molecule.
        remove_parent (bool, optional): Remove parent from the fragments.
        sanitize (bool, optional): Wether to sanitize the fragments.
        fix (bool, optional): Wether to fix the fragments.
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
    mol: Chem.Mol,
    remove_parent: bool = False,
    sanitize: bool = True,
    fix: bool = True,
):
    """Fragment molecule by applying brics first, then fall back to frag.

    Args:
        mols (Chem.Mol): a molecule.
        remove_parent (bool, optional): Remove parent from the fragments.
        sanitize (bool, optional): Wether to sanitize the fragments.
        fix (bool, optional): Wether to fix the fragments.
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
    mol,
    pattern: str = None,
    max_cut: int = 1,
    max_bond_cut: int = 20,
    h_split: bool = False,
):
    """Fragment molecule on specific bonds suitable for a MMPA analysis.

    Args:
        mol (Chem.Mol): Molecule to fragment.
        pattern (str, optional): Bond pattern to split on.
            Will use default rdkit pattern '[#6+0;!$(*=,#[!#6])]!@!=!#[*]' if not provided
        max_cut (int, optional): Number of cuts. Default to 3.
        max_bond_cut (int, optional): Maximum number of bond to cut. Default to 20.
        h_split (bool, optional):  Whether to split at hydrogen position too.
            This is equivalent to enabling the addition of new fragments.
            Default to False.

    Returns:
        List of fragments
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
        mol = Chem.AddHs(mol)
        frags += rdMMPA.FragmentMol(
            mol,
            pattern="[#1]!@!=!#[!#1]",
            maxCuts=1,
            resultsAsMols=False,
            maxCutBonds=max_bond_cut,
        )
    return set(frags)


def mmpa_cut(mol: Chem.Mol, rdkit_pattern: bool = False):
    """Cut molecules to perform mmpa analysis later

    Args:
        mol (Chem.Mol or str): Molecule to fragment.
        rdkit_pattern (bool, optional): Whether to perform the fragmentation
            using the default rdkit pattern: [#6+0;!$(*=, #[!#6])]!@!=!#[*]"
            Default to False.

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
    mol = Chem.AddHs(mol)
    smiles = dm.to_smiles(mol)

    n = mol.GetNumHeavyAtoms()
    if n < 60:
        frags = mmpa_frag(mol, pattern=None, max_cut=1, max_bond_cut=100, h_split=True)
        for core, chains in frags:
            output = f"{smiles},{core},{chains}\n"
            outlines.add(output)

    return outlines
