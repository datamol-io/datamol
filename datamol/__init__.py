from typing import TYPE_CHECKING

import os

import importlib


# The below lazy import logic is coming from openff-toolkit:
# https://github.com/openforcefield/openff-toolkit/blob/b52879569a0344878c40248ceb3bd0f90348076a/openff/toolkit/__init__.py#L44

# Dictionary of objects to lazily import; maps the object's name to its module path
_lazy_imports_obj = {
    # version
    "__version__": "datamol._version",
    "is_lower_than_current_rdkit_version": "datamol._version",
    "is_greater_than_current_rdkit_version": "datamol._version",
    "is_lower_eq_than_current_rdkit_version": "datamol._version",
    "is_greater_eq_than_current_rdkit_version": "datamol._version",
    # types
    "Mol": "datamol.types",
    "BondType": "datamol.types",
    "ChemicalReaction": "datamol.types",
    "Atom": "datamol.types",
    "Bond": "datamol.types",
    "DatamolColor": "datamol.types",
    "RDKitColor": "datamol.types",
    # utils
    "parallelized": "datamol.utils",
    "parallelized_with_batches": "datamol.utils",
    "JobRunner": "datamol.utils",
    "fs": "datamol.utils",
    # log
    "freesolv": "datamol.data",
    "cdk2": "datamol.data",
    "solubility": "datamol.data",
    "chembl_drugs": "datamol.data",
    "chembl_samples": "datamol.data",
    # log
    "enable_rdkit_log": "datamol.log",
    "disable_rdkit_log": "datamol.log",
    "without_rdkit_log": "datamol.log",
    # mol
    "PERIODIC_TABLE": "datamol.mol",
    "TRIPLE_BOND": "datamol.mol",
    "DOUBLE_BOND": "datamol.mol",
    "SINGLE_BOND": "datamol.mol",
    "AROMATIC_BOND": "datamol.mol",
    "UNSPECIFIED_BOND": "datamol.mol",
    "copy_mol": "datamol.mol",
    "to_mol": "datamol.mol",
    "same_mol": "datamol.mol",
    "reorder_atoms": "datamol.mol",
    "randomize_atoms": "datamol.mol",
    "to_neutral": "datamol.mol",
    "sanitize_mol": "datamol.mol",
    "sanitize_first": "datamol.mol",
    "sanitize_smiles": "datamol.mol",
    "standardize_smiles": "datamol.mol",
    "standardize_mol": "datamol.mol",
    "fix_valence_charge": "datamol.mol",
    "incorrect_valence": "datamol.mol",
    "decrease_bond": "datamol.mol",
    "fix_valence": "datamol.mol",
    "adjust_singleton": "datamol.mol",
    "remove_dummies": "datamol.mol",
    "fix_mol": "datamol.mol",
    "replace_dummies_atoms": "datamol.mol",
    "keep_largest_fragment": "datamol.mol",
    "is_transition_metal": "datamol.mol",
    "set_dative_bonds": "datamol.mol",
    "set_mol_props": "datamol.mol",
    "copy_mol_props": "datamol.mol",
    "atom_indices_to_mol": "datamol.mol",
    "protect_atoms": "datamol.mol",
    "atom_list_to_bond": "datamol.mol",
    "substructure_matching_bonds": "datamol.mol",
    "add_hs": "datamol.mol",
    "remove_hs": "datamol.mol",
    "unique_id": "datamol.mol",
    "hash_mol": "datamol.mol",
    "clear_mol_props": "datamol.mol",
    "strip_mol_to_core": "datamol.mol",
    "make_scaffold_generic": "datamol.mol",
    "to_scaffold_murcko": "datamol.mol",
    "compute_ring_system": "datamol.mol",
    "clear_atom_props": "datamol.mol",
    "clear_atom_map_number": "datamol.mol",
    "set_atom_positions": "datamol.mol",
    "get_atom_positions": "datamol.mol",
    "remove_salts_solvents": "datamol.mol",
    # cluster
    "cluster_mols": "datamol.cluster",
    "pick_diverse": "datamol.cluster",
    "pick_centroids": "datamol.cluster",
    "assign_to_centroids": "datamol.cluster",
    # convert
    "to_smiles": "datamol.convert",
    "to_selfies": "datamol.convert",
    "from_selfies": "datamol.convert",
    "to_smarts": "datamol.convert",
    "from_smarts": "datamol.convert",
    "smiles_as_smarts": "datamol.convert",
    "to_inchi": "datamol.convert",
    "to_inchikey": "datamol.convert",
    "from_inchi": "datamol.convert",
    "to_df": "datamol.convert",
    "from_df": "datamol.convert",
    "render_mol_df": "datamol.convert",
    "to_inchi_non_standard": "datamol.convert",
    "to_inchikey_non_standard": "datamol.convert",
    # fp
    "to_fp": "datamol.fp",
    "fp_to_array": "datamol.fp",
    "list_supported_fingerprints": "datamol.fp",
    "fold_count_fp": "datamol.fp",
    # similarity
    "pdist": "datamol.similarity",
    "cdist": "datamol.similarity",
    # io
    "read_csv": "datamol.io",
    "read_excel": "datamol.io",
    "read_sdf": "datamol.io",
    "to_sdf": "datamol.io",
    "to_smi": "datamol.io",
    "read_smi": "datamol.io",
    "read_mol2file": "datamol.io",
    "read_molblock": "datamol.io",
    "to_molblock": "datamol.io",
    "to_xlsx": "datamol.io",
    "read_pdbblock": "datamol.io",
    "to_pdbblock": "datamol.io",
    "read_pdbfile": "datamol.io",
    "to_pdbfile": "datamol.io",
    "open_df": "datamol.io",
    "save_df": "datamol.io",
    # isomers
    "enumerate_stereoisomers": "datamol.isomers",
    "enumerate_tautomers": "datamol.isomers",
    "enumerate_structisomers": "datamol.isomers",
    "canonical_tautomer": "datamol.isomers",
    "remove_stereochemistry": "datamol.isomers",
    # viz
    "to_image": "datamol.viz",
    "lasso_highlight_image": "datamol.viz",
    # mcs
    "find_mcs": "datamol.mcs",
    # graph
    "to_graph": "datamol.graph",
    "get_all_path_between": "datamol.graph",
    "match_molecular_graphs": "datamol.graph",
    "reorder_mol_from_template": "datamol.graph",
    # Remember to add new lazy imports to __all__ and the if TYPE_CHECKING imports
}

# Dictionary of modules to lazily import; maps the modules's name to its path
_lazy_imports_mod = {
    "fragment": "datamol.fragment",
    "scaffold": "datamol.scaffold",
    "molar": "datamol.molar",
    "descriptors": "datamol.descriptors",
    "predictors": "datamol.predictors",
    "reactions": "datamol.reactions",
    "convert": "datamol.convert",
    "fp": "datamol.fp",
    "similarity": "datamol.similarity",
    "io": "datamol.io",
    "isomers": "datamol.isomers",
    "mcs": "datamol.mcs",
    "graph": "datamol.graph",
    "align": "datamol.align",
    "viz": "datamol.viz",
    "conformers": "datamol.conformers",
    "utils": "datamol.utils",
    "data": "datamol.data",
}


def __getattr__(name):
    """Lazily import objects from _lazy_imports_obj or _lazy_imports_mod

    Note that this method is only called by Python if the name cannot be found
    in the current module."""
    obj_mod = _lazy_imports_obj.get(name)
    if obj_mod is not None:
        mod = importlib.import_module(obj_mod)
        return mod.__dict__[name]

    lazy_mod = _lazy_imports_mod.get(name)
    if lazy_mod is not None:
        return importlib.import_module(lazy_mod)

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    """Add _lazy_imports_obj and _lazy_imports_mod to dir(<module>)"""
    keys = (*globals().keys(), *_lazy_imports_obj.keys(), *_lazy_imports_mod.keys())
    return sorted(keys)


if TYPE_CHECKING or os.environ.get("DATAMOL_DISABLE_LAZY_LOADING", "0") == "1":
    # These types are imported lazily at runtime, but we need to tell type
    # checkers what they are.

    from ._version import __version__
    from ._version import is_lower_than_current_rdkit_version
    from ._version import is_greater_than_current_rdkit_version
    from ._version import is_lower_eq_than_current_rdkit_version
    from ._version import is_greater_eq_than_current_rdkit_version

    from .types import Mol
    from .types import BondType
    from .types import ChemicalReaction
    from .types import Atom
    from .types import Bond
    from .types import DatamolColor
    from .types import RDKitColor

    from . import utils

    from .utils import parallelized
    from .utils import parallelized_with_batches
    from .utils import JobRunner
    from .utils import fs

    from .data import freesolv
    from .data import cdk2
    from .data import solubility
    from .data import chembl_drugs
    from .data import chembl_samples

    from .log import enable_rdkit_log
    from .log import disable_rdkit_log
    from .log import without_rdkit_log

    from .mol import PERIODIC_TABLE
    from .mol import TRIPLE_BOND
    from .mol import DOUBLE_BOND
    from .mol import SINGLE_BOND
    from .mol import AROMATIC_BOND
    from .mol import UNSPECIFIED_BOND

    from .mol import copy_mol
    from .mol import to_mol
    from .mol import same_mol
    from .mol import reorder_atoms
    from .mol import randomize_atoms
    from .mol import to_neutral
    from .mol import sanitize_mol
    from .mol import sanitize_first
    from .mol import sanitize_smiles
    from .mol import standardize_smiles
    from .mol import standardize_mol
    from .mol import fix_valence_charge
    from .mol import incorrect_valence
    from .mol import decrease_bond
    from .mol import fix_valence
    from .mol import adjust_singleton
    from .mol import remove_dummies
    from .mol import fix_mol
    from .mol import replace_dummies_atoms
    from .mol import keep_largest_fragment
    from .mol import is_transition_metal
    from .mol import set_dative_bonds
    from .mol import set_mol_props
    from .mol import copy_mol_props
    from .mol import atom_indices_to_mol
    from .mol import protect_atoms
    from .mol import atom_list_to_bond
    from .mol import substructure_matching_bonds
    from .mol import add_hs
    from .mol import remove_hs
    from .mol import unique_id
    from .mol import hash_mol
    from .mol import clear_mol_props
    from .mol import strip_mol_to_core
    from .mol import make_scaffold_generic
    from .mol import to_scaffold_murcko
    from .mol import compute_ring_system
    from .mol import clear_atom_props
    from .mol import clear_atom_map_number
    from .mol import set_atom_positions
    from .mol import get_atom_positions
    from .mol import remove_salts_solvents

    from .cluster import cluster_mols
    from .cluster import pick_diverse
    from .cluster import pick_centroids
    from .cluster import assign_to_centroids

    from .convert import to_smiles
    from .convert import to_selfies
    from .convert import from_selfies
    from .convert import to_smarts
    from .convert import from_smarts
    from .convert import smiles_as_smarts
    from .convert import to_inchi
    from .convert import to_inchikey
    from .convert import from_inchi
    from .convert import to_df
    from .convert import from_df
    from .convert import render_mol_df
    from .convert import to_inchi_non_standard
    from .convert import to_inchikey_non_standard

    from .fp import to_fp
    from .fp import fp_to_array
    from .fp import list_supported_fingerprints
    from .fp import fold_count_fp

    from .similarity import pdist
    from .similarity import cdist

    from .io import read_csv
    from .io import read_excel
    from .io import read_sdf
    from .io import to_sdf
    from .io import to_smi
    from .io import read_smi
    from .io import read_mol2file
    from .io import read_molblock
    from .io import to_molblock
    from .io import to_xlsx
    from .io import read_pdbblock
    from .io import to_pdbblock
    from .io import read_pdbfile
    from .io import to_pdbfile
    from .io import save_df
    from .io import open_df

    from .isomers import enumerate_stereoisomers
    from .isomers import enumerate_tautomers
    from .isomers import enumerate_structisomers
    from .isomers import canonical_tautomer
    from .isomers import remove_stereochemistry

    from . import align
    from . import conformers
    from . import viz
    from . import fragment
    from . import scaffold
    from . import molar
    from . import descriptors
    from . import predictors
    from . import reactions

    from .viz import to_image
    from .viz import lasso_highlight_image

    from .mcs import find_mcs

    from .graph import to_graph
    from .graph import get_all_path_between
    from .graph import match_molecular_graphs
    from .graph import reorder_mol_from_template
