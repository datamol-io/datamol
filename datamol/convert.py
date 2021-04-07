from typing import Union
from typing import List
from typing import Optional

import re

from loguru import logger

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
        cxsmiles: Whether to return a CXSMILES instead of a SMILES.
        allow_to_fail: Raise an error if the conversion to SMILES fails. Return None otherwise.
    """
    if ordered and canonical is False:
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

    selfies = sf.encoder(mol)  # type: ignore

    if selfies == -1:
        return None

    return selfies


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
    inchi: Optional[str],
    sanitize: bool = True,
    remove_hs: bool = True,
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
    render_df_mol: bool = True,
    render_all_df_mol: bool = False,
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
        render_df_mol: whether to render the molecule in the dataframe to images.
            If called once, it will be applied for the newly created dataframe with
            mol in it.
        render_all_df_mol: Whether to render all pandas dataframe mol column as images.
    """

    # Init a dataframe
    df = pd.DataFrame()

    # Feed it with smiles
    if smiles_column is not None:
        smiles = [dm.to_smiles(mol) for mol in mols]
        df[smiles_column] = smiles

    # Add a mol column
    if mol_column is not None:
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

    if smiles_column is not None and smiles_column in props_df.columns:
        logger.warning(
            f"The SMILES column name provided ('{smiles_column}') is already present in the properties"
            " of the molecules. THe returned dataframe will two columns with the same name."
        )

    # Concat the df with the properties df
    df = pd.concat([df, props_df], axis=1)

    # Render mol column to images
    if render_df_mol is True and mol_column is not None:
        # NOTE(hadim): replace by `PandaTools.ChangeMoleculeRendering` once
        # https://github.com/rdkit/rdkit/issues/3563 is fixed.
        _ChangeMoleculeRendering(df)

        if render_all_df_mol:
            PandasTools.RenderImagesInAllDataFrames()

    return df


def from_df(
    df: pd.DataFrame,
    smiles_column: Optional[str] = "smiles",
    mol_column: str = None,
    conserve_smiles: bool = False,
) -> List[Chem.rdchem.Mol]:
    """Convert a dataframe to a list of mols.

    Note:
        If `smiles_column` is used to build the molecules, this property
        is removed from the molecules' properties. You can decide to conserve
        the SMILES column by setting `conserve_smiles` to True.

    Args:
        df: a dataframe.
        smiles_column: Column name to extract the molecule.
        mol_column: Column name to extract the molecule. It takes
            precedence over `smiles_column`.
        conserve_smiles: Whether to conserve the SMILES in the mols' props.
    """

    if smiles_column is None and mol_column is None:
        raise ValueError("Either `smiles_column` or `mol_column` must be not None.")

    if len(df) == 0:
        return []

    def _row_to_mol(row):

        props = row.to_dict()

        if mol_column is not None:
            mol = props[mol_column]
        else:

            if conserve_smiles:
                smiles = props[smiles_column]
            else:
                # If a SMILES column is used to create the molecule then it is removed from the
                # properties.
                smiles = props.pop(smiles_column)

            mol = dm.to_mol(smiles)

        if mol is None:
            return None

        dm.set_mol_props(mol, props)
        return mol

    return df.apply(_row_to_mol, axis=1).tolist()


def _ChangeMoleculeRendering(frame=None, renderer="PNG"):
    """Allows to change the rendering of the molecules between base64 PNG images and string
    representations.
    This serves two purposes: First it allows to avoid the generation of images if this is
    not desired and, secondly, it allows to enable image rendering for newly created dataframe
    that already contains molecules, without having to rerun the time-consuming
    AddMoleculeColumnToFrame. Note: this behaviour is, because some pandas methods, e.g. head()
    returns a new dataframe instance that uses the default pandas rendering (thus not drawing
    images for molecules) instead of the monkey-patched one.
    """
    import types

    if renderer == "String":
        Chem.rdchem.Mol.__str__ = PandasTools.PrintDefaultMolRep
    else:
        Chem.rdchem.Mol.__str__ = PandasTools.PrintAsBase64PNGString

    if frame is not None:
        frame.to_html = types.MethodType(PandasTools.patchPandasHTMLrepr, frame)

    if PandasTools.defPandasRepr is not None and renderer == "String":
        frame._repr_html_ = types.MethodType(PandasTools.defPandasRepr, frame)
    else:
        frame._repr_html_ = types.MethodType(PandasTools.patchPandasrepr, frame)
