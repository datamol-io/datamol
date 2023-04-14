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

from . import utils

from .utils import parallelized
from .utils import parallelized_with_batches
from .utils import JobRunner
from .utils import fs

from .data import freesolv
from .data import cdk2
from .data import solubility

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

from .cluster import cluster_mols
from .cluster import pick_diverse
from .cluster import pick_centroids
from .cluster import assign_to_centroids

from . import fragment
from . import scaffold
from . import molar
from . import descriptors
from . import predictors
from . import reactions

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

from .graph import to_graph
from .graph import get_all_path_between

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

from .isomers import enumerate_stereoisomers
from .isomers import enumerate_tautomers
from .isomers import enumerate_structisomers
from .isomers import canonical_tautomer
from .isomers import remove_stereochemistry

from . import align
from . import conformers

# from . import viz

# from .viz import to_image
# from .viz import lasso_highlight_image

from .mcs import find_mcs

from .graph import to_graph
from .graph import get_all_path_between
from .graph import match_molecular_graphs
from .graph import reorder_mol_from_template

import importlib

# Dictionary of objects to lazily import; maps the object's name to its module path
_lazy_imports_obj = {
    "to_image": "datamol.viz",
    "lasso_highlight_image": "datamol.viz",
    # Remember to add new lazy imports to __all__ and the if TYPE_CHECKING imports
}

# Dictionary of modules to lazily import; maps the modules's name to its path
_lazy_imports_mod = {
    "viz": "datamol.viz",
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
