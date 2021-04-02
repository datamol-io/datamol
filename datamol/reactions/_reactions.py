from functools import singledispatch

import numpy as np
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem

import datamol as dm


def is_reaction_ok(rxn):
    """Check is the reaction is synthetically valid"""
    return rdChemReactions.SanitizeRxn(rxn) == rdChemReactions.SanitizeFlags.SANITIZE_NONE


def compute_reaction_product(out, single_output=True):
    """Compute the product of a reaction"""
    out = [dm.fix_mol(x[0], n_iter=0) for x in out]
    if not single_output:
        return [dm.sanitize_mol(x) for x in out]
    # Might be a important to make a tradeoff decision in selecting products for greater speed.
    # product = sorted(out, key=lambda x: MoleculeEnv.compute_reward_from_mol(x, True))[-1]
    # sampling from list of products is an alternative
    return dm.sanitize_first(np.random.permutation(out))


@singledispatch
def apply_reaction(rxn, mol, react_pos, single_output=False):
    """Apply a chemical reaction on a molecule"""
    raise ValueError


@apply_reaction.register(rdChemReactions.ChemicalReaction)
def _(rxn, mol, react_pos, single_output=False):
    # only cares about the first product right now
    # Anyway, there is only one major product in the database

    if not rxn.IsInitialized():
        rxn.Initialize()
    out = rxn.RunReactant(mol, react_pos)
    return compute_reaction_product(out, single_output)


@apply_reaction.register(tuple)
def _(rxn, mol, react_pos, single_output=False):
    reaction, reactants = rxn
    if not reaction.IsInitialized():
        reaction.Initialize()
    # now substitute one of the reactant by the mol
    reactants = list(reactants)
    reactants[react_pos] = mol
    out = reaction.RunReactants(tuple(reactants))
    return compute_reaction_product(out, single_output)


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
