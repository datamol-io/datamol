# `datamol`

Datamol is designed to be used with a single import (`import datamol as dm`). Most of the functions are available in `datamol.*`. The others ones are available throught their specific modules.

The below sections shows you the directly available Datamol functions. For other modules, please browser the API using the left menu.

## Working with molecules

### The basics

::: datamol:to_mol
::: datamol:copy_mol
::: datamol:reorder_atoms
::: datamol:randomize_atoms
::: datamol:to_neutral
::: datamol:set_mol_props
::: datamol:copy_mol_props
::: datamol:atom_indices_to_mol

### Fix, sanitize and standardize

::: datamol:sanitize_mol
::: datamol:sanitize_first
::: datamol:sanitize_smiles
::: datamol:standardize_smiles
::: datamol:standardize_mol
::: datamol:fix_valence_charge
::: datamol:incorrect_valence
::: datamol:decrease_bond
::: datamol:fix_valence
::: datamol:adjust_singleton
::: datamol:remove_dummies
::: datamol:fix_mol
::: datamol:replace_dummies_atoms
::: datamol:keep_largest_fragment
::: datamol:is_transition_metal
::: datamol:set_dative_bonds

### Enumerate

::: datamol:enumerate_stereoisomers
::: datamol:enumerate_tautomers

## Convert molecule(s)

::: datamol:to_smiles
::: datamol:to_selfies
::: datamol:from_selfies
::: datamol:to_smarts
::: datamol:to_inchi
::: datamol:to_inchikey
::: datamol:from_inchi
::: datamol:to_df
::: datamol:from_df

## Input/Output

::: datamol:read_csv
::: datamol:read_excel
::: datamol:read_sdf
::: datamol:to_sdf
::: datamol:to_smi
::: datamol:read_smi

## Molecule similarity and distance

::: datamol:pdist
::: datamol:cdist

## Working with fingerprints

::: datamol:to_fp
::: datamol:fp_to_array
::: datamol:list_supported_fingerprints
::: datamol:fold_count_fp

## Cluster molecules

::: datamol:cluster_mols
::: datamol:pick_diverse
::: datamol:pick_centroids
::: datamol:assign_to_centroids

## Molecule as a graph

::: datamol:to_graph
::: datamol:get_all_path_between

## Constants

::: datamol:PERIODIC_TABLE
::: datamol:TRIPLE_BOND
::: datamol:DOUBLE_BOND
::: datamol:SINGLE_BOND
::: datamol:AROMATIC_BOND

## Control RDKit logging

::: datamol:without_rdkit_log
::: datamol:enable_rdkit_log
::: datamol:disable_rdkit_log

## Toy dataset

::: datamol:freesolv
