import pytest

import datamol as dm
import numpy as np


def test_esol():
    smiles = "Nc1cnn(-c2ccccc2)c(=O)c1Cl"
    mol = dm.to_mol(smiles)

    assert np.allclose(dm.predictors.esol(mol), -2.627091966265316)


def test_esol_from_data():
    data = dm.freesolv()
    data = data.iloc[:20]

    with pytest.raises(KeyError):
        dm.predictors.esol_from_data(data)

    data["mol"] = data["smiles"].apply(dm.to_mol)
    data["clogp"] = data["mol"].apply(dm.descriptors.clogp)
    data["mw"] = data["mol"].apply(dm.descriptors.mw)
    data["n_rotatable_bonds"] = data["mol"].apply(dm.descriptors.n_rotatable_bonds)
    data["n_aromatic_atoms_proportion"] = data["mol"].apply(
        dm.descriptors.n_aromatic_atoms_proportion
    )

    # dataframe
    esol_values = dm.predictors.esol_from_data(data)
    assert esol_values.dtype == float
    assert esol_values.shape == (20,)

    # series
    v = dm.predictors.esol_from_data(data.iloc[0])
    v = float(v)
    assert isinstance(v, float)

    # dict
    v = dm.predictors.esol_from_data(data.iloc[0].to_dict())
    v = float(v)
    assert isinstance(v, float)
