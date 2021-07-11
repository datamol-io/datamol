from typing import Union
from typing import Optional

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdmolops
from rdkit.Chem import rdReducedGraphs
from rdkit.Chem.EState import Fingerprinter as EStateFingerprinter
from rdkit.Avalon import pyAvalonTools

from rdkit.DataStructs.cDataStructs import ConvertToNumpyArray
from rdkit.DataStructs.cDataStructs import SparseBitVect
from rdkit.DataStructs.cDataStructs import UIntSparseIntVect
from rdkit.DataStructs.cDataStructs import IntSparseIntVect
from rdkit.DataStructs.cDataStructs import LongSparseIntVect
from rdkit.DataStructs.cDataStructs import ExplicitBitVect
from rdkit.DataStructs.cDataStructs import ULongSparseIntVect

import numpy as np

import datamol as dm

_FP_FUNCS = {
    "maccs": rdMolDescriptors.GetMACCSKeysFingerprint,
    "ecfp": rdMolDescriptors.GetMorganFingerprintAsBitVect,
    "topological": rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect,
    "atompair": rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect,
    "rdkit": rdmolops.RDKFingerprint,
    "pattern": rdmolops.PatternFingerprint,
    "layered": rdmolops.LayeredFingerprint,
    "erg": rdReducedGraphs.GetErGFingerprint,
    # NOTE(hadim): bad for pickling?
    "estate": lambda x, **args: EStateFingerprinter.FingerprintMol(x)[0],
    "avalon-count": pyAvalonTools.GetAvalonCountFP,
    # NOTE(hadim): `rdkit-count` to ndarray hangs forever on my machine. Disabling it for now.
    # "rdkit-count": rdmolops.UnfoldedRDKFingerprintCountBased,
    "ecfp-count": rdMolDescriptors.GetHashedMorganFingerprint,
    "fcfp-count": rdMolDescriptors.GetHashedMorganFingerprint,
    "topological-count": rdMolDescriptors.GetHashedTopologicalTorsionFingerprint,
    "atompair-count": rdMolDescriptors.GetHashedAtomPairFingerprint,
}


def fp_to_array(
    fp: Union[np.ndarray, SparseBitVect, ExplicitBitVect, UIntSparseIntVect]
) -> np.ndarray:
    """Convert rdkit fingerprint to numpy array.

    Note:
        This implementation has shown to be faster than using `DataStructs.ConvertToNumpyArray`
        by a factor of ~4. See https://github.com/rdkit/rdkit/discussions/3863.

    Args:
        fp: The fingerprint.
    """

    if isinstance(fp, np.ndarray):
        fp_out = fp

    elif isinstance(fp, SparseBitVect):
        tmp = np.zeros(fp.GetNumBits(), dtype=int)
        for n_bit in list(fp.GetOnBits()):
            tmp[n_bit] = 1
        fp_out = tmp

    elif isinstance(fp, ExplicitBitVect):
        fp_out = np.frombuffer(fp.ToBitString().encode(), "u1") - ord("0")

    elif isinstance(
        fp, (UIntSparseIntVect, IntSparseIntVect, LongSparseIntVect, ULongSparseIntVect)
    ):
        # one of the other rdkit type
        tmp = np.array(0, dtype="int")
        ConvertToNumpyArray(fp, tmp)
        fp_out = tmp

    else:
        raise ValueError(
            f"The fingerprint of type '{type(fp)}' is not supported. "
            "Please open a ticket at https://github.com/datamol-org/datamol/issues."
        )

    return fp_out


def to_fp(
    mol: Union[str, Chem.rdchem.Mol],
    fp_size: int = 2048,
    radius: int = 3,
    use_features: bool = True,
    fp_args: dict = None,
    as_array: bool = True,
    fp_type: str = "ecfp",
) -> Optional[Union[np.ndarray, SparseBitVect, ExplicitBitVect, UIntSparseIntVect]]:
    """Transform a molecule from smiles to morgan fingerprint.

    Note:
        That function should be expanded to compute more type of fingerprints.

    Args:
        mol: a molecule or a SMILES.
        fp_size: Size of morgan fingerprint. Default to 2048.
            Warning: only apply to the default ECFP fingerprint.
        radius: Radius of the morgan fingerprints. Default to 3.
            Warning: only apply to the default ECFP fingerprint.
        use_features: Whether to use atom features. Default to True.
            Warning: only apply to the default ECFP fingerprint.
        fp_args: Arguments to build the fingerprint. Refer to the official RDKit documentation.
        as_array: Whether to return a numpy array of an RDKit vec. Default to True.
        fp_type: The type of the fingerprint. See `dm.list_supported_fingerprints()` for a
            complete list.

    Returns:
        A fingerprint vector or None
    """

    # Get fp function
    fp_func = _FP_FUNCS.get(fp_type)

    if fp_func is None:
        raise ValueError(
            f"The fingerprint '{fp_type}' is not available. Use `dm.list_supported_fingerprints()` to "
            "get a complete list of the available fingerprints."
        )

    # Convert input to mol if needed
    if isinstance(mol, str):
        mol = dm.to_mol(mol)

    if fp_args is None:
        fp_args = {}

    # Use the default fp args only for the "ECFP" fingerprint.
    # NOTE(hadim): not sure what to do here in term of backward compatibility. The current
    # solution keeps it but leads to an inconsistent API.
    # NOTE(hadim): also should I use `FP_DEF_PARAMS`? Or let the default being the RDKit ones? I think
    # I prefer the latter one but some args are required (`radius` for example).
    if fp_type in ["ecfp", "ecfp-count", "fcfp-count"]:
        if "fp_size" not in fp_args:
            fp_args["nBits"] = fp_size
        if "radius" not in fp_args:
            fp_args["radius"] = radius
        if "use_features" not in fp_args:
            fp_args["useFeatures"] = use_features

    fp = fp_func(mol, **fp_args)

    if as_array:
        return fp_to_array(fp)

    return fp


def list_supported_fingerprints():
    """Return the supported fingerprints in datamol."""

    return set(_FP_FUNCS.keys())
