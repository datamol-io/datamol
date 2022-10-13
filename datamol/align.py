from typing import Optional
from typing import Union
from typing import Sequence
from typing import Any

from packaging import version
from collections import defaultdict as ddict

import rdkit
from rdkit.Chem import rdDepictor
from rdkit.Chem import rdMolAlign

import pandas as pd
import datamol as dm


from .types import Mol


def compute_2d_coords(mol: Mol, copy: bool = True, verbose: bool = False) -> Mol:
    """Compute 2D coordinates for a molecule.

    Args:
        mol: A molecule.
        copy: Whether to copy the molecule.
    """
    if copy:
        mol = dm.copy_mol(mol)

    with dm.without_rdkit_log(enable=not verbose):
        rdDepictor.Compute2DCoords(mol)

    return mol


def template_align(
    mol: Union[str, Mol],
    template: Optional[Union[str, Mol]] = None,
    copy: bool = True,
    use_depiction: bool = True,
    remove_confs: bool = True,
    auto_select_coord_gen: bool = False,
) -> Optional[Mol]:
    """Align an input molecule to a template. If the template is not provided then the input molecule is
    returned.

    Args:
        mol: A molecule.
        template: Template to align to.
        copy: whether to copy the molecule before aligning it.
        use_depiction: Whether to use the depiction API or use MolAlign
            The main difference is around how 3D information is handled, but also, because the depiction API
            will emphasize the atoms that do not match, whereas AlignMol will not.
        remove_confs: Whether to remove all conformation in the input molecule first.
            You can set this to true when not using depiction
        auto_select_coord_gen: Whether to automatically select the coordinate generation method.

    Returns:
        mol: aligned molecule (dm.Mol). None if initial mol argument is undefined or invalid.
    """

    if isinstance(mol, str):
        _mol = dm.to_mol(mol)
    elif copy:
        _mol = dm.copy_mol(mol)
    else:
        _mol = mol

    if _mol is None:
        return None

    if isinstance(template, str):
        _template = dm.to_mol(template)
    elif copy:
        _template = dm.copy_mol(template)
    else:
        _template = template

    if _template is None:
        return _mol

    if remove_confs:
        _mol.RemoveAllConformers()

    # EN: to make this more general and robust, we need to first check whether the template
    # has 2D coordinates or not. If it does not, we need to compute them.
    # This is a very rare edge case if requests are comming from the UI, but better to check than not xD
    if _template.GetNumConformers() == 0:
        _template = compute_2d_coords(_template)

    # EN: now we can align the molecule to the template
    # but first, we should avoid MCS as much as possible, because it's expensive
    # so if the template is a subgraph of the molecule, no need to perform any MCS
    # Another reason for this, is to avoid inconsistency in alignment between molecules that are
    # supergraph of the template vs molecules that are subgraph of the template.
    pattern = _template
    if not _mol.HasSubstructMatch(_template):

        pattern = None
        mcs_smarts = dm.find_mcs([_mol, _template])

        if mcs_smarts is not None:
            pattern = dm.from_smarts(mcs_smarts)

    if pattern is not None:

        if auto_select_coord_gen:
            rdDepictor.SetPreferCoordGen(use_depiction)

        # we would need to compute 2d coordinates for the molecules if it doesn't have any
        if _mol.GetNumConformers() == 0:
            _mol = compute_2d_coords(_mol)

        if use_depiction:
            rdDepictor.GenerateDepictionMatching2DStructure(
                _mol,
                reference=_template,
                refPatt=pattern,
                acceptFailure=True,
                allowRGroups=True,
            )
        else:
            query_match = _mol.GetSubstructMatch(pattern)
            template_match = _template.GetSubstructMatch(pattern)
            rdMolAlign.AlignMol(_mol, _template, atomMap=list(zip(query_match, template_match)))

    return _mol


def auto_align_many(
    mols: Union[Sequence[Mol], pd.Series],
    partition_method: str = "anon-scaffold",
    copy: bool = True,
    cluster_cutoff: float = 0.7,
    allow_r_groups: bool = True,
    **kwargs: Any,
):
    """Partition a list of molecules into clusters sharing common scaffold of common core,
    then align the molecules to that common core. This function will compute the list of
    smiles/smarts representative of each cluster first.

    The returned molecules will have two properties associated to them:

    - `dm.auto_align_many.cluster_id`: the cluster id of the molecule.
    - `dm.auto_align_many.core`: the smiles/smarts of the core of the cluster.

    Args:
        mols: A list of molecules to auto align.
        partition_method: Partition method to use:

            - 'scaffold': Cluster molecules by Murcko scaffold.
            - 'strip-scaffold': Cluster molecules by Murcko scaffold, but remove all atoms not
                in the core.
            - 'anon-scaffold': Cluster molecules by Murcko scaffold, but making it
                generic including the bonds.
            - 'anongraph-scaffold': Cluster molecules by Murcko scaffold, but making it
                generic but keeping the bond order informations.
            - 'cluster': Cluster the molecules using Butina frm RDKit with `dm.cluster_mols`.
            Cautious as the method 'cluster' is very sensitive to the cutoff.

        copy: Whether to copy the molecules before aligning them.
        cluster_cutoff: Optional cluster cutoff.
        allow_r_groups: Optional, if True, terminal dummy atoms in the
                        reference are ignored if they match an implicit hydrogen in the
                        molecule, and a constrained depiction is still attempted
        **kwargs: Additional arguments to pass to clustering method
    """

    if copy:
        mols = [dm.copy_mol(mol) for mol in mols]

    mol_groups = ddict(list)  # map scaffold index to list of unique molecules
    scaffold_mols = {}

    if partition_method.endswith("scaffold"):

        scaffolds = [dm.to_scaffold_murcko(m) for m in mols]
        scaffolds_ids = [dm.to_smiles(x) for x in scaffolds]

        if partition_method.startswith("strip-"):
            scaffolds = [dm.strip_mol_to_core(x) for x in scaffolds]
            scaffolds_ids = [dm.to_smiles(x) for x in scaffolds]

        elif partition_method.startswith("anongraph-"):
            scaffolds = [dm.make_scaffold_generic(s, include_bonds=True) for s in scaffolds]
            scaffolds_ids = [dm.to_smiles(x) for x in scaffolds]

        elif partition_method.startswith("anon-"):
            scaffolds = [dm.make_scaffold_generic(s, include_bonds=False) for s in scaffolds]
            scaffolds_ids = [dm.to_smiles(x) for x in scaffolds]

        for i, s in enumerate(scaffolds_ids):
            mol_groups[s].append(i)
            scaffolds[i] = compute_2d_coords(scaffolds[i])
            scaffold_mols[s] = scaffolds[i]

    elif partition_method == "cluster":
        # partition is cluster, first compute molecule clusters
        clusters, mol_clusters = dm.cluster_mols(mols, cutoff=cluster_cutoff, **kwargs)

        # now compute the mcs for each clusters
        cluster_mcs = [
            (dm.find_mcs(mol_cluster) if len(mol_cluster) > 1 else dm.to_smiles(mol_cluster[0]))
            for mol_cluster in mol_clusters
        ]
        scaffolds_ids = [cluster_mcs[cluster_id] for cluster_id, _ in enumerate(clusters)]

        for i, s in enumerate(scaffolds_ids):
            mol_groups[s].extend(clusters[i])

        for x in scaffolds_ids:
            core = None
            if x is not None:
                core = dm.from_smarts(x)
                core = compute_2d_coords(core)
            scaffold_mols[x] = core

    else:
        raise ValueError(f"Unknown partition method: {partition_method}")

    # now we match each molecule to the scaffold and align them
    # note that the molecule object will be modified in place in the list
    for cluster_id, (core, mols_ids) in enumerate(mol_groups.items()):
        core_mol = scaffold_mols[core]

        for mol_id in mols_ids:

            mol = mols[mol_id]

            if core_mol is not None:
                # Only pass allowRGroups if current
                # rdkit version is >= 2021_03_1
                # ref https://github.com/rdkit/rdkit/pull/3811
                allowRGroups = (
                    {"allowRGroups": allow_r_groups}
                    if version.parse(rdkit.__version__) >= version.parse("2021.03.1")
                    else {}
                )
                rdDepictor.GenerateDepictionMatching2DStructure(
                    mol,
                    reference=core_mol,
                    acceptFailure=True,
                    **allowRGroups,
                )

            # Add some props to the mol so the user can retrieve the groups from
            # it later.
            props = {}
            props["dm.auto_align_many.cluster_id"] = cluster_id
            props["dm.auto_align_many.core"] = core
            dm.set_mol_props(mol, props)

    # EN: you can discard the mol_groups (or keep it and match the values
    # to molecular line notation, so you will not have to reocompute the above)
    # and convert the mols into cxsmiles if you want
    # return mols, mol_groups

    return mols
