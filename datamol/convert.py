from typing import Union
from typing import List
from typing import Optional

import re

from loguru import logger

import pandas as pd

from rdkit import Chem
from rdkit.Chem import rdmolfiles
from rdkit.Chem import PandasTools

import selfies as sf

from .types import Mol

# NOTE(hadim): it's not possible to use `from .mol import ...` because of circular imports
import datamol as dm


def to_smiles(
    mol: Mol,
    canonical: bool = True,
    isomeric: bool = True,
    kekulize: bool = False,
    ordered: bool = False,
    explicit_bonds: bool = False,
    explicit_hs: bool = False,
    randomize: bool = False,
    cxsmiles: bool = False,
    allow_to_fail: bool = False,
    with_atom_indices: bool = False,
) -> Optional[str]:
    """Convert a mol to a SMILES.

    Args:
        mol: a molecule.
        canonical: if false no attempt will be made to canonicalize the molecule.
        isomeric: whether to include information about stereochemistry in the SMILES.
        kekulize: whether to return the kekule version of the SMILES.
        ordered: whether to force reordering of the atoms first.
        explicit_bonds: if true, all bond orders will be explicitly indicated in the output SMILES.
        explicit_hs: if true, all H counts will be explicitly indicated in the output SMILES.
        randomize: whether to randomize the generated smiles. Override `canonical`.
        cxsmiles: Whether to return a CXSMILES instead of a SMILES.
        allow_to_fail: Raise an error if the conversion to SMILES fails. Return None otherwise.
        with_atom_indices: Whether to add atom indices to the SMILES.
    """

    if ordered and canonical is False:
        mol = dm.reorder_atoms(mol)

    if randomize:
        mol = dm.randomize_atoms(mol)
        canonical = False

    if with_atom_indices:
        mol = dm.atom_indices_to_mol(mol, copy=True)

    smiles = None
    try:

        if cxsmiles:
            smiles = rdmolfiles.MolToCXSmiles(
                mol,
                isomericSmiles=isomeric,
                canonical=canonical,
                allBondsExplicit=explicit_bonds,
                allHsExplicit=explicit_hs,
                kekuleSmiles=kekulize,
            )

        else:
            smiles = rdmolfiles.MolToSmiles(
                mol,
                isomericSmiles=isomeric,
                canonical=canonical,
                allBondsExplicit=explicit_bonds,
                allHsExplicit=explicit_hs,
                kekuleSmiles=kekulize,
            )

    except Exception as e:

        if allow_to_fail:
            raise e

        return None

    return smiles


def to_selfies(mol: Union[str, Mol]) -> Optional[str]:
    """Convert a mol to SELFIES.

    Args:
        mol: a molecule or a SMILES.

    Returns:
        selfies: SELFIES string.
    """

    if isinstance(mol, Mol):
        mol = to_smiles(mol)

    if mol is None:
        return None
    selfies = sf.encoder(mol)

    if selfies == -1:
        return None

    return selfies


def from_selfies(selfies: str, as_mol: bool = False) -> Optional[Union[str, Mol]]:
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


def smiles_as_smarts(mol: Union[str, Mol], keep_hs: bool = True) -> Optional[str]:
    """Convert a smiles to a smarts if possible

    Args:
        mol: a molecule.
        keep_hs: Whether to keep hydrogen. This will increase the count of H atoms
            for atoms with attached hydrogens to create a valid smarts without further substitution allowed
            e.g. [H]-[CH]-[*] -> [H]-[CH2]-[*]

    Returns:
        smarts of the molecule
    """

    if isinstance(mol, str):
        mol = dm.to_mol(mol)

    if mol is None:
        return None

    # Change the isotope to 99
    for atom in mol.GetAtoms():  # type: ignore
        if keep_hs:
            s = sum(na.GetAtomicNum() == 1 for na in atom.GetNeighbors())
            if s:
                atom.SetNumExplicitHs(atom.GetTotalNumHs() + s)
        atom.SetIsotope(99)

    # Print out the smiles, all the atom attributes will be fully specified
    smarts = to_smiles(mol, isomeric=True, explicit_bonds=True)

    if smarts is None:
        return None

    # Remove the 99 isotope labels
    smarts = re.sub(r"\[99", "[", smarts)
    return smarts


def to_inchi(mol: Union[str, Mol]) -> Optional[str]:
    """Convert a mol to a standard Inchi.

    Args:
        mol: a molecule.
    """

    if isinstance(mol, str):
        mol = dm.to_mol(mol)

    if mol is None:
        return None

    inchi_val = Chem.MolToInchi(mol)
    if not inchi_val:
        return None
    return inchi_val


def to_inchi_non_standard(
    mol: Union[str, Mol],
    fixed_hydrogen_layer: bool = True,
    undefined_stereocenter: bool = True,
    reconnected_metal_layer: bool = True,
    tautomerism_keto_enol: bool = True,
    tautomerism_15: bool = True,
    options: Optional[List[str]] = None,
) -> Optional[str]:
    """Convert a mol to a non-standard Inchi.

    Note that turning all the flags to `False` will result in the standard Inchi.

    **Warning**: this function will return a **non-standard** Inchi. See
    https://www.inchi-trust.org/technical-faq-2 for details.

    It's important to not mix standard and non-standard InChi. If you don't know
    much about non-standard InChi, we highly recommend you to use the
    standard InChi with `dm.to_inchi()`.

    Args:
        mol: a molecule.
        fixed_hydrogen_layer: whether to include a fixed hydrogen layer (`/FixedH`).
        undefined_stereocenter: whether to include an undefined stereocenter layer (`/SUU`).
        reconnected_metal_layer: whether to include reconnected metals (`/RecMet`).
        tautomerism_keto_enol: whether to account tautomerism keto-enol (`/KET`).
        tautomerism_15: whether to account 1,5-tautomerism (`/15T`).
        options: More InchI options in a form of a list of string. Example:
            `["/SRel", "/AuxNone"]`.
    """

    if isinstance(mol, str):
        mol = dm.to_mol(mol)

    if mol is None:
        return None

    inchi_options = _process_inchi_options(
        fixed_hydrogen_layer=fixed_hydrogen_layer,
        undefined_stereocenter=undefined_stereocenter,
        reconnected_metal_layer=reconnected_metal_layer,
        tautomerism_keto_enol=tautomerism_keto_enol,
        tautomerism_15=tautomerism_15,
        options=options,
    )

    inchi_val = Chem.MolToInchi(mol, options=inchi_options)
    if not inchi_val:
        return None
    return inchi_val


def to_smarts(mol: Mol) -> Optional[str]:
    """Convert a mol to SMARTS format

    Args:
        mol: a molecule.
    """

    if mol is None:
        return None

    return Chem.MolToSmarts(mol)  # type: ignore


def to_inchikey(mol: Union[str, Mol]) -> Optional[str]:
    """Convert a mol to a standard InchiKey.

    Args:
        mol: a molecule
    """
    if isinstance(mol, str):
        mol = dm.to_mol(mol)

    if mol is None:
        return None

    inchikey = Chem.MolToInchiKey(mol)
    if not inchikey:
        return None
    return inchikey


def to_inchikey_non_standard(
    mol: Union[str, Mol],
    fixed_hydrogen_layer: bool = True,
    undefined_stereocenter: bool = True,
    reconnected_metal_layer: bool = True,
    tautomerism_keto_enol: bool = True,
    tautomerism_15: bool = True,
    options: Optional[List[str]] = None,
) -> Optional[str]:
    """Convert a mol to a non-standard InchiKey.

    Note that turning all the flags to `False` will result in the standard InchiKey.

    **Warning**: this function will return a **non-standard** InchiKey. See
    https://www.inchi-trust.org/technical-faq-2 for details.

    It's important to not mix standard and non-standard InChiKey. If you don't know
    much about non-standard InchiKey, we highly recommend you to use the
    standard InchiKey with `dm.to_inchikey()`.

    Args:
        mol: a molecule
        fixed_hydrogen_layer: whether to include a fixed hydrogen layer (`/FixedH`).
        undefined_stereocenter: whether to include an undefined stereocenter layer (`/SUU`).
        reconnected_metal_layer: whether to include reconnected metals (`/RecMet`).
        tautomerism_keto_enol: whether to account tautomerism keto-enol (`/KET`).
        tautomerism_15: whether to account 1,5-tautomerism (`/15T`).
        options: More InchI options in a form of a list of string. Example:
            `["/SRel", "/AuxNone"]`.
    """

    if isinstance(mol, str):
        mol = dm.to_mol(mol)

    if mol is None:
        return None

    inchi_options = _process_inchi_options(
        fixed_hydrogen_layer=fixed_hydrogen_layer,
        undefined_stereocenter=undefined_stereocenter,
        reconnected_metal_layer=reconnected_metal_layer,
        tautomerism_keto_enol=tautomerism_keto_enol,
        tautomerism_15=tautomerism_15,
        options=options,
    )

    inchikey = Chem.MolToInchiKey(mol, options=inchi_options)
    if not inchikey:
        return None
    return inchikey


def from_inchi(
    inchi: Optional[str],
    sanitize: bool = True,
    remove_hs: bool = True,
) -> Optional[Mol]:
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


def from_smarts(smarts: Optional[str]) -> Optional[Mol]:
    """Convert a SMARTS string to a molecule

    Args:
        smarts: a smarts string
    """

    if smarts is None:
        return None
    return Chem.MolFromSmarts(smarts)  # type: ignore


def to_df(
    mols: List[Mol],
    smiles_column: Optional[str] = "smiles",
    mol_column: Optional[str] = None,
    include_private: bool = False,
    include_computed: bool = False,
    render_df_mol: bool = True,
    render_all_df_mol: bool = False,
) -> Optional[pd.DataFrame]:
    """Convert a list of mols to a dataframe using each mol properties
    as a column.

    For the reverse operation, you might to check `dm.from_df()`.

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
        smiles = [to_smiles(mol) for mol in mols]
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

        render_mol_df(df)

        if render_all_df_mol:
            PandasTools.RenderImagesInAllDataFrames()

    return df


def from_df(
    df: pd.DataFrame,
    smiles_column: Optional[str] = "smiles",
    mol_column: Optional[str] = None,
    conserve_smiles: bool = False,
    sanitize: bool = True,
) -> List[Mol]:
    """Convert a dataframe to a list of mols.

    For the reverse operation, you might to check `dm.to_df()`.

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
        sanitize: Whether to sanitize if `smiles_column` is not None.
    """

    if smiles_column is None and mol_column is None:
        raise ValueError("Either `smiles_column` or `mol_column` must be not None.")

    if len(df) == 0:
        return []

    # Try to detect the mol column if `mol_column` is None.
    if mol_column is None:
        for col in df.columns:
            if isinstance(df[col].iloc[0], Mol):
                mol_column = col

    def _row_to_mol(row):

        props = row.to_dict()

        if mol_column is not None:
            mol = props.pop(mol_column)
        else:

            if conserve_smiles:
                smiles = props[smiles_column]
            else:
                # If a SMILES column is used to create the molecule then it is removed from the
                # properties.
                smiles = props.pop(smiles_column)

            mol = dm.to_mol(smiles, sanitize=sanitize)

        if mol is None:
            return None

        dm.set_mol_props(mol, props)
        return mol

    return df.apply(_row_to_mol, axis=1).tolist()


def render_mol_df(df: pd.DataFrame):
    """Render the molecules column in a dataframe. The rendering is performed
    in-place only. So nothing is returned.

    Args:
        df: a dataframe.
    """
    # NOTE(hadim): replace by `PandaTools.ChangeMoleculeRendering` once
    # https://github.com/rdkit/rdkit/issues/3563 is fixed.
    _ChangeMoleculeRendering(df)


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
        frame._repr_html_ = types.MethodType(PandasTools.defPandasRepr, frame)  # type: ignore
    else:
        frame._repr_html_ = types.MethodType(PandasTools.patchPandasrepr, frame)  # type: ignore


def _process_inchi_options(
    fixed_hydrogen_layer: bool = True,
    undefined_stereocenter: bool = True,
    reconnected_metal_layer: bool = True,
    tautomerism_keto_enol: bool = True,
    tautomerism_15: bool = True,
    options: Optional[List[str]] = None,
):

    inchi_options = []

    if fixed_hydrogen_layer:
        inchi_options.append("/FixedH")

    if undefined_stereocenter:
        inchi_options.append("/SUU")

    if reconnected_metal_layer:
        inchi_options.append("/RecMet")

    if tautomerism_keto_enol:
        inchi_options.append("/KET")

    if tautomerism_15:
        inchi_options.append("/15T")

    if options is not None:
        inchi_options.extend(options)

    inchi_options = " ".join(inchi_options)
    return inchi_options
