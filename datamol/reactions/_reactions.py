from typing import cast
from typing import Union

import os
import io

from loguru import logger
import numpy as np

import fsspec

from rdkit.Chem import rdChemReactions

import datamol as dm


ATTACHING_RXN = rdChemReactions.ReactionFromSmarts("[*;h:1]>>[*:1][*]")


def rxn_from_smarts(rxn_smarts: str) -> dm.ChemicalReaction:
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


def rxn_to_smarts(rxn: dm.ChemicalReaction) -> str:
    """Create a SMARTS from a reaction.

    Args:
        rxn: dm.ChemicalReaction object.

    Returns:
        SMARTS as string.
    """
    return rdChemReactions.ReactionToSmarts(reaction=rxn)


def rxn_from_block(
    rxn_block: str,
    sanitize: bool = False,
) -> dm.ChemicalReaction:
    """Create a reaction from a block.

    Args:
        rxn_block: A reaction block.
        sanitize: Whether to sanitize the reaction.

    Returns:
        Initialized reaction.

    """
    rxn = rdChemReactions.ReactionFromRxnBlock(rxnblock=rxn_block, sanitize=sanitize)
    rxn.Initialize()
    return rxn


def rxn_from_block_file(
    rxn_block_path: Union[str, os.PathLike],
    sanitize: bool = False,
) -> dm.ChemicalReaction:
    """Create a reaction from a block file.

    Args:
        rxn_block_path: Filepath to a reaction block file.
        sanitize: Whether to sanitize the reaction.

    Returns:
        Initialized reaction.
    """
    with fsspec.open(rxn_block_path) as f:
        rxn_block = f.read()  # type: ignore
        rxn = rxn_from_block(rxn_block=rxn_block, sanitize=sanitize)
    return rxn


def rxn_to_block(
    rxn: dm.ChemicalReaction,
    separate_agents: bool = False,
    force_V3000: bool = False,
) -> str:
    """Create a block from a reaction object.

    Args:
        rxn: A reaction object.
        separate_agents: Whether to separate agents from the reactants block. Not supported
            if `force_V3000=True`.
        force_V3000: Write the block in a V3000 format.

    Returns:
        Reaction block as string
    """

    args = {}
    if dm.is_lower_than_current_rdkit_version("2022"):
        logger.warning("RDKit version prior to 2022.* does not support the `force_V3000` flag.")
    else:
        args["forceV3000"] = force_V3000

    return rdChemReactions.ReactionToRxnBlock(reaction=rxn, separateAgents=separate_agents, **args)


def rxn_to_block_file(
    rxn: dm.ChemicalReaction,
    output_block_path: Union[str, os.PathLike],
    separate_agents: bool = False,
    force_V3000: bool = False,
):
    """Create a block from a reaction object.

    Args:
        rxn: A reaction object.
        output_block_path: Filepath to a reaction block file.
        separate_agents: Whether to separate agents from the reactants block. Not supported
            if `force_V3000=True`.
        force_V3000: Write the block in a V3000 format.
    """
    block = rxn_to_block(
        rxn=rxn,
        separate_agents=separate_agents,
        force_V3000=force_V3000,
    )

    with fsspec.open(output_block_path, "w") as f:
        f = cast(io.TextIOWrapper, f)
        f.write(block)


def is_reaction_ok(rxn: dm.ChemicalReaction, enable_logs: bool = False) -> bool:
    """Check if the given reaction is synthetically valid.

    Args:
        rxn: dm.ChemicalReaction object
        enable_logs: Whether to enable logs.

    Returns:
        Boolean whether reaction is valid
    """
    nWarn, nError, nReactants, nProducts, labels = rdChemReactions.PreprocessReaction(rxn)

    if enable_logs:
        logger.info(f"Number of warnings:{nWarn}")
        logger.info(f"Number of preprocessing errors: {nError}")
        logger.info(f"Number of reactants in reaction: {nReactants}")
        logger.info(f"Number of products in reaction: {nProducts}")
        logger.info(f"Preprocess labels added:{labels}")

    return rdChemReactions.SanitizeRxn(rxn) in [
        rdChemReactions.SanitizeFlags.SANITIZE_NONE,
        None,
    ]


def select_reaction_output(
    product: list,
    single_output: bool = True,
    rm_attach: bool = False,
    as_smiles: bool = False,
    sanitize: bool = True,
) -> Union[list, str, dm.Mol]:
    """
    Compute the products from a reaction. It only takes the first product of the

    Args:
        product: All the products from a reaction.
        single_output: Whether return a single output from a reaction.
        rm_attach: Whether remove the attachment point from the product.
        as_smiles: Whether return the result in smiles.
        sanitize: Whether sanitize the product to return.

    Returns:
        Processed products from reaction.
    """
    # flatten all possible products of a reaction
    product = list(sum(product, ()))
    if single_output:
        product = list(np.random.choice(product, 1))
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
    rxn: dm.ChemicalReaction,
    reactants: tuple,
    single_output: bool = False,
    as_smiles: bool = False,
    rm_attach: bool = False,
    disable_logs: bool = True,
    sanitize: bool = True,
) -> Union[list, str, dm.Mol]:
    """
    Apply a chemical reaction on a molecule

    Args:
       rxn: Reaction object.
       reactants: A tuple of reactants.
       single_output: Whether return one product from all possible product.
       as_smiles: Whether return product in SMILES.
       rm_attach: Whether remove the attachment point from product.
       disable_logs: Whether disable rdkit logs.
       sanitize: Whether sanitize the product.

    Returns:
       Reaction products.
    """
    with dm.without_rdkit_log(enable=disable_logs):
        if not rxn.IsInitialized():
            rxn.Initialize()  # pragma: no cover

        product = rxn.RunReactants(reactants)
        outputs = select_reaction_output(
            product=product,
            single_output=single_output,
            as_smiles=as_smiles,
            rm_attach=rm_attach,
            sanitize=sanitize,
        )

    return outputs


def can_react(rxn: dm.ChemicalReaction, mol: dm.Mol) -> bool:
    """Check if a molecule is a reactant to a chemical reaction.

    Args:
        rxn: Reaction to check.
        mol: Molecule to check if it is a reactant.

    Returns:
        True if `mol` is a reactant of rxn.
    """
    if not rxn.IsInitialized():
        rxn.Initialize()  # pragma: no cover
    if rxn.IsMoleculeReactant(mol):
        return find_reactant_position(rxn, mol) != -1
    return False


def find_reactant_position(rxn: dm.ChemicalReaction, mol: dm.Mol) -> int:
    """Find the position of a reactant in a reaction.

    Args:
        rxn: Reaction
        mol: Molecule

    Returns:
        Reactant position or -1 if `mol` is not a reactant.
    """

    if not rxn.IsInitialized():
        rxn.Initialize()  # pragma: no cover

    react_pos = -1
    for pos, rct in enumerate(rxn.GetReactants()):
        if mol.HasSubstructMatch(rct):
            react_pos = pos
    return react_pos


def inverse_reaction(rxn: dm.ChemicalReaction) -> dm.ChemicalReaction:
    """
    Get the reverse reaction of the input reaction

    Args:
        rxn: Reaction to inverse.

    Returns:
        Inversed reaction.
    """
    rxn2 = rdChemReactions.ChemicalReaction()
    for i in range(rxn.GetNumReactantTemplates()):
        rxn2.AddProductTemplate(rxn.GetReactantTemplate(i))
    for i in range(rxn.GetNumProductTemplates()):
        rxn2.AddReactantTemplate(rxn.GetProductTemplate(i))
    rxn2.Initialize()
    return rxn2
