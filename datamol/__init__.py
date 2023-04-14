import importlib


from ._version import __version__

from ._version import is_lower_than_current_rdkit_version
from ._version import is_greater_than_current_rdkit_version
from ._version import is_lower_eq_than_current_rdkit_version
from ._version import is_greater_eq_than_current_rdkit_version


# Dictionary of objects to lazily import; maps the object's name to its module path
_lazy_imports_obj = {
    # types
    "Mol": "datamol.types",
    "BondType": "datamol.types",
    "ChemicalReaction": "datamol.types",
    "Atom": "datamol.types",
    "Bond": "datamol.types",
    # utils
    "parallelized": "datamol.utils",
    "parallelized_with_batches": "datamol.utils",
    "JobRunner": "datamol.utils",
    "fs": "datamol.utils",
    # log
    "freesolv": "datamol.data",
    "cdk2": "datamol.data",
    "solubility": "datamol.data",
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
