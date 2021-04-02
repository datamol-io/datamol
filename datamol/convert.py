from typing import Union
from typing import List
from typing import Optional

import re

import pandas as pd

from rdkit import Chem
from rdkit.Chem import PandasTools

import selfies as sf

import datamol as dm


def to_smiles(
    mol: Chem.rdchem.Mol,
    canonical: bool = True,
    isomeric: bool = True,
    ordered: bool = False,
    explicit_bonds: bool = False,
    explicit_hs: bool = False,
    randomize: bool = False,
    cxsmiles: bool = False,
    allow_to_fail: bool = False,
) -> Optional[str]:
    """Convert a mol to a SMILES.

    Args:
        mol: a molecule.
        canonical: if false no attempt will be made to canonicalize the molecule.
        isomeric: whether to include information about stereochemistry in the SMILES.
        ordered: whether to force reordering of the atoms first.
        explicit_bonds: if true, all bond orders will be explicitly indicated in the output SMILES.
        explicit_hs: if true, all H counts will be explicitly indicated in the output SMILES.
        randomize: whether to randomize the generated smiles. Override `canonical`.
    """
    if ordered:
        mol = dm.reorder_atoms(mol)

    if randomize:
        mol = dm.randomize_atoms(mol)
        canonical = False

    smiles = None
    try:

        if cxsmiles:
            smiles = Chem.MolToCXSmiles(  # type: ignore
                mol,
                isomericSmiles=isomeric,
                canonical=canonical,
                allBondsExplicit=explicit_bonds,
                allHsExplicit=explicit_hs,
            )

        else:
            smiles = Chem.MolToSmiles(  # type: ignore
                mol,
                isomericSmiles=isomeric,
                canonical=canonical,
                allBondsExplicit=explicit_bonds,
                allHsExplicit=explicit_hs,
            )

    except Exception as e:

        if allow_to_fail:
            raise e

        return None

    return smiles


def to_selfies(mol: Union[str, Chem.rdchem.Mol]) -> Optional[str]:
    """Convert a mol to SELFIES.

    Args:
        mol: a molecule or a SMILES.

    Returns:
        selfies: SELFIES string.
    """
    if mol is None:
        return None

    if isinstance(mol, Chem.rdchem.Mol):
        mol = to_smiles(mol)

    return sf.encoder(mol)


def from_selfies(selfies: str, as_mol: bool = False) -> Optional[Union[str, Chem.rdchem.Mol]]:
    """Convert a SEFLIES to a smiles or a mol.

    Args:
        selfies: a selfies.
        as_mol (str, optional): whether to return a mol or a smiles.

    Returns:
        smiles or mol.
    """
    if selfies is None:
        return None

    smiles = sf.decoder(selfies)

    if as_mol and smiles is not None:
        return dm.to_mol(smiles)

    return smiles


def to_smarts(mol: Union[str, Chem.rdchem.Mol], keep_hs: bool = True) -> Optional[str]:
    """Convert a molecule to a smarts.

    Args:
        mol: a molecule.
        keep_hs: Whether to keep hydrogen. This will increase the count of H atoms
            for atoms with attached hydrogens to create a valid smarts.
            e.g. [H]-[CH2]-[*] -> [H]-[CH3]-[*]

    Returns:
        smarts of the molecule
    """

    if mol is None:
        return None

    if isinstance(mol, str):
        mol = dm.to_mol(mol)

    # Change the isotope to 42
    for atom in mol.GetAtoms():  # type: ignore
        if keep_hs:
            s = sum(na.GetAtomicNum() == 1 for na in atom.GetNeighbors())
            if s:
                atom.SetNumExplicitHs(atom.GetTotalNumHs() + s)
        atom.SetIsotope(42)

    # Print out the smiles, all the atom attributes will be fully specified
    smarts = to_smiles(mol, isomeric=True, explicit_bonds=True)

    if smarts is None:
        return None

    # Remove the 42 isotope labels
    smarts = re.sub(r"\[42", "[", smarts)
    return smarts


def to_inchi(mol: Union[str, Chem.rdchem.Mol]) -> Optional[str]:
    """Convert a mol to Inchi.

    Args:
        mol: a molecule.
    """

    if mol is None:
        return None

    if isinstance(mol, str):
        mol = dm.to_mol(mol)

    return Chem.MolToInchi(mol)


def to_inchikey(mol: Union[str, Chem.rdchem.Mol]) -> Optional[str]:
    """Convert a mol to Inchi key.

    Args:
        mol: a molecule
    """

    if mol is None:
        return None

    if isinstance(mol, str):
        mol = dm.to_mol(mol)

    return Chem.MolToInchiKey(mol)


def from_inchi(
    inchi: str, sanitize: bool = True, remove_hs: bool = True
) -> Optional[Chem.rdchem.Mol]:
    """Convert an InChi to a mol.

    Args:
        inchi: an inchi string.
        sanitize: do sanitize.
        remove_hs: do remove hs.

    Returns:
        mol
    """
    if inchi is None:
        return None

    return Chem.MolFromInchi(inchi, sanitize=sanitize, removeHs=remove_hs)


def to_df(
    mols: List[Chem.rdchem.Mol],
    smiles_column: Optional[str] = "smiles",
    mol_column: str = None,
    include_private: bool = False,
    include_computed: bool = False,
) -> Optional[pd.DataFrame]:
    """Convert a list of mols to a dataframe using each mol properties
    as a column.

    Args:
        mols: a molecule.
        smiles_column: name of the SMILES column.
        mol_column: Name of the column. If not None, rdkit.Chem.PandaTools
            is used to add a molecule column.
        include_private: Include private properties in the columns.
        include_computed: Include computed properties in the columns.
    """

    # Init a dataframe
    df = pd.DataFrame()

    # Feed it with smiles
    if smiles_column is not None:
        smiles = [dm.to_smiles(mol) for mol in mols]
        df[smiles_column] = smiles

    # Add a mol column
    if mol_column is not None:
        # NOTE(hadim): not sure it's ok to do that here.
        PandasTools.RenderImagesInAllDataFrames()
        df[mol_column] = mols

    # Add any other properties present in the molecule
    props = [
        mol.GetPropsAsDict(
            includePrivate=include_private,
            includeComputed=include_computed,
        )
        for mol in mols
    ]
    props_df = pd.DataFrame(props)

    # If a smiles column already exists in the props, we remove it
    # in favor of the one specified as args.
    if smiles_column is not None and smiles_column in props_df.columns:
        props_df.pop(smiles_column)

    # Concat the df with the properties df
    df = pd.concat([df, props_df], axis=1)

    return df


def from_df(
    df: pd.DataFrame,
    smiles_column: Optional[str] = "smiles",
    mol_column: str = None,
) -> List[Chem.rdchem.Mol]:
    """Convert a dataframe to a list of mols.

    Args:
        df: a dataframe.
        smiles_column: Column name to extract the molecule.
        mol_column: Column name to extract the molecule. It takes
            precedence over `smiles_column`.
    """

    if smiles_column is None and mol_column is None:
        raise ValueError("Either `smiles_column` or `mol_column` must be not None.")

    if len(df) == 0:
        return []

    def _row_to_mol(row):

        if mol_column is not None:
            mol = row[mol_column]
        else:
            mol = dm.to_mol(row[smiles_column])

        if mol is None:
            return None

        dm.set_mol_props(mol, row.to_dict())
        return mol

    return df.apply(_row_to_mol, axis=1).tolist()
