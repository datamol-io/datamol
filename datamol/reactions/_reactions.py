from functools import singledispatch
from typing import Union
from loguru import logger
import numpy as np

import fsspec

from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem

import datamol as dm


ATTACHING_RXN = rdChemReactions.ReactionFromSmarts("[*;h:1]>>[*:1][*]")


def rxn_from_smarts(rxn_smarts: str) -> rdChemReactions.ChemicalReaction:
    """
    Create a reaction from smarts

    Args:
        rxn_smarts:  Reaction SMARTS string

    Returns:
        Initilized reaction.
    """
    rxn = rdChemReactions.ReactionFromSmarts(SMARTS=rxn_smarts)
    rxn.Initialize()
    return rxn


def rxn_from_block(
    rxn_block: str = None, rxn_file: str = None, sanitize=False
) -> rdChemReactions.ChemicalReaction:
    """
    Create a reaction from smarts

    Args:
        rxn_block_or_file:  Reaction SMARTS string or file path

    Returns:
        Initilized reaction.
    """
    if rxn_file is not None and rxn_block is None:
        with fsspec.open(rxn_block_or_file) as f:
            rxn_block = f.read(f)

    rxn = rdChemReactions.ReactionFromRxnBlock(rxnblock=rxn_block, sanitize=sanitize)
    rxn.Initialize()
    return rxn


def is_reaction_ok(rxn: Union[str, rdChemReactions.ChemicalReaction]) -> bool:
    """
    Check if the given reaction is synthetically valid

    Args:
        rxn: rdkit.rdChemReactions.ChemicalReaction object

    Returns:
        Boolean whether reaction is valid
    """
    nWarn, nError, nReactants, nProducts, labels = AllChem.PreprocessReaction(rxn)
    logger.info(f"Number of warnings:{nWarn}")
    logger.info(f"Number of preprocessing errors: {nError}")
    logger.info(f"Number of reactants in reaction: {nReactants}")
    logger.info(f"Number of products in reaction: {nProducts}")
    logger.info(f"Preprocess labels added:{labels}")
    return rdChemReactions.SanitizeRxn(rxn) == rdChemReactions.SanitizeFlags.SANITIZE_NONE


def compute_reaction_product(
    product,
    single_output=True,
    rm_attach: bool = False,
    as_smiles: bool = False,
    sanitize: bool = True,
) -> Union[list, str, Chem.Mol]:
    """
    Compute the products from a reaction.

    Args:
        product: All the products from a reaction.
        single_output: Whether return a single output from a reaction.
        rm_attach: Whether remove the attachment point from the product.
        as_smiles: Whether return the result in smiles.
        sanitize: Whether sanitize the product to return.

    Returns:
        Processed products from reaction.
    """
    product = [x[0] for x in product]
    if single_output:
        product = np.random.choice(product, 1)
    if sanitize:
        product = [dm.sanitize_mol(m) for m in product]
    if rm_attach:
        product = [dm.remove_dummies(x) for x in product]
    if as_smiles:
        product = [dm.to_smiles(x) for x in product if x is not None]
    if single_output:
        return product[0]
    return product


def apply_reaction(
    rxn: rdChemReactions.ChemicalReaction,
    reactants: tuple,
    single_output: bool = False,
    as_smiles: bool = False,
    rm_attach: bool = False,
    disable_log: bool = True,
    sanitize: bool = True,
) -> Union[list, str, Chem.Mol]:
    """
    Apply a chemical reaction on a molecule

    Args:
       rxn: Reaction object.
       reactants: A tuple of reactants.
       single_output: Whether return one product from all possible product.
       as_smiles: Whether return product in SMILES.
       rm_attach: Whether remove the attachment point from product.
       disable_log: Whether disable rdkit log.
       sanitize: Whether sanitize the product.

    Returns:
       Reaction products.
    """
    if disable_log:
        dm.disable_rdkit_log()
    if not rxn.IsInitialized():
        rxn.Initialize()
    product = rxn.RunReactants(reactants)
    return compute_reaction_product(
        product=product,
        single_output=single_output,
        as_smiles=as_smiles,
        rm_attach=rm_attach,
        sanitize=sanitize,
    )


def can_react(rxn, mol):
    """
    Check if a molecule is a reactant to a chemical reaction and return position

    Args:
        rxn: Reaction to check.
        mol: Molecule to check if it is a reactant.

    Returns:
        Position of mol
    """
    if not rxn.IsInitialized():
        rxn.Initialize()
    if rxn.IsMoleculeReactant(mol):
        return _find_rct_position(rxn, mol)
    return False


def _find_rct_position(rxn, mol):
    """
    Find the position of a reactant in a reaction

    Args:
        rxn: Reaction
        mol: Molecule

    Returns:
        Reactant position
    """
    react_pos = -1
    for pos, rct in enumerate(rxn.GetReactants()):
        if mol.HasSubstructMatch(rct):
            react_pos = pos
    return react_pos


def inverse_reaction(rxn):
    """
    Get the reverse reaction of the input reaction

    Args:
        rxn: Reaction to inverse.

    Returns:
        Inversed reaction.
    """
    rxn2 = AllChem.ChemicalReaction()
    for i in range(rxn.GetNumReactantTemplates()):
        rxn2.AddProductTemplate(rxn.GetReactantTemplate(i))
    for i in range(rxn.GetNumProductTemplates()):
        rxn2.AddReactantTemplate(rxn.GetProductTemplate(i))
    rxn2.Initialize()
    return rxn2
