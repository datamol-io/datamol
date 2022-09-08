from functools import singledispatch
from typing import Union
import numpy as np
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem

import datamol as dm


ATTACHING_RXN = rdChemReactions.ReactionFromSmarts("[*;h:1]>>[*:1][*]")


def from_smarts(rxn_smarts: str):
    """Create a reaction from smarts"""
    return rdChemReactions.ReactionFromSmarts(rxn_smarts)


def is_reaction_ok(rxn: Union[str, rdChemReactions.ChemicalReaction]):
    """Check is the reaction is synthetically valid"""
    if isinstance(rxn, str):
        rxn = from_smarts(rxn_smarts=rxn)
    return rdChemReactions.SanitizeRxn(rxn) == rdChemReactions.SanitizeFlags.SANITIZE_NONE


def compute_reaction_product(
    out, single_output=True, rm_attach: bool = False, as_smiles: bool = False, sanitize: bool = True
):
    """Compute the products of a reaction"""
    out = [x[0] for x in out]
    if single_output:
        out = np.random.choice(out, 1)
    if sanitize:
        out = [dm.sanitize_mol(m) for m in out]
    if rm_attach:
        out = [dm.remove_dummies(x) for x in out]
    if as_smiles:
        out = [dm.to_smiles(x) for x in out if x is not None]
    if single_output:
        return out[0]
    return out

    # Might be a important to make a tradeoff decision in selecting products for greater speed.
    # product = sorted(out, key=lambda x: MoleculeEnv.compute_reward_from_mol(x, True))[-1]
    # sampling from list of products is an alternative


def apply_reaction(
    rxn: rdChemReactions.ChemicalReaction,
    reactants: tuple,
    single_output: bool = False,
    as_smiles: bool = False,
    rm_attach: bool = False,
    disable_log: bool = True,
    sanitize: bool = True,
):
    """Apply a chemical reaction on a molecule"""
    if disable_log:
        dm.disable_rdkit_log()
    if not rxn.IsInitialized():
        rxn.Initialize()
    out = rxn.RunReactants(reactants)
    return compute_reaction_product(
        out=out,
        single_output=single_output,
        as_smiles=as_smiles,
        rm_attach=rm_attach,
        sanitize=sanitize,
    )


@singledispatch
def can_react(rxn, mol):
    """Check if a molecule is a reactant to a chemical reaction and return position"""
    raise ValueError


@can_react.register(rdChemReactions.ChemicalReaction)
def _(rxn, mol):
    if not rxn.IsInitialized():
        rxn.Initialize()
    return _find_rct_position(rxn, mol)


@can_react.register(tuple)
def _(rxn, mol):
    reaction = rxn[0]
    if not reaction.IsInitialized():
        reaction.Initialize()
    return _find_rct_position(reaction, mol)


def _find_rct_position(rxn, mol):
    """Find the position of a reactant in a reaction"""
    react_pos = -1
    for pos, rct in enumerate(rxn.GetReactants()):
        if mol.HasSubstructMatch(rct):
            react_pos = pos
    return react_pos


def inverse_reaction(rxn):
    """Get the reverse reaction of the input reaction"""
    rxn2 = AllChem.ChemicalReaction()
    for i in range(rxn.GetNumReactantTemplates()):
        rxn2.AddProductTemplate(rxn.GetReactantTemplate(i))
    for i in range(rxn.GetNumProductTemplates()):
        rxn2.AddReactantTemplate(rxn.GetProductTemplate(i))
    rxn2.Initialize()
    return rxn2
