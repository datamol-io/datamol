import pytest

import pandas as pd
import datamol as dm
import numpy as np


def test_descriptors():

    mol = dm.to_mol("CCN(CC)CCCC(C)NC1=C2C=CC(=CC2=NC=C1)Cl")
    dm.descriptors.mw(mol)
    dm.descriptors.fsp3(mol)
    dm.descriptors.n_hba(mol)
    dm.descriptors.n_hbd(mol)
    dm.descriptors.n_lipinski_hba(mol)
    dm.descriptors.n_lipinski_hbd(mol)
    dm.descriptors.n_rings(mol)
    dm.descriptors.n_hetero_atoms(mol)
    dm.descriptors.n_heavy_atoms(mol)
    dm.descriptors.n_rotatable_bonds(mol)
    dm.descriptors.n_aliphatic_rings(mol)
    dm.descriptors.n_aromatic_rings(mol)
    dm.descriptors.n_saturated_rings(mol)
    dm.descriptors.n_radical_electrons(mol)
    dm.descriptors.tpsa(mol)
    dm.descriptors.qed(mol)
    dm.descriptors.clogp(mol)
    dm.descriptors.sas(mol)


def test_compute_many_descriptors():

    mol = dm.to_mol("CCN(CC)CCCC(C)NC1=C2C=CC(=CC2=NC=C1)Cl")

    true_values = pd.Series(
        {
            "mw": 319.181525512,
            "fsp3": 0.5,
            "n_lipinski_hba": 3.0,
            "n_lipinski_hbd": 1.0,
            "n_rings": 2.0,
            "n_hetero_atoms": 4.0,
            "n_heavy_atoms": 22.0,
            "n_rotatable_bonds": 8.0,
            "n_radical_electrons": 0.0,
            "tpsa": 28.16,
            "qed": 0.7564117572128701,
            "clogp": 4.810600000000004,
            "sas": 2.670786229594949,
            "n_aliphatic_carbocycles": 0.0,
            "n_aliphatic_heterocyles": 0.0,
            "n_aliphatic_rings": 0.0,
            "n_aromatic_carbocycles": 1.0,
            "n_aromatic_heterocyles": 1.0,
            "n_aromatic_rings": 2.0,
            "n_saturated_carbocycles": 0.0,
            "n_saturated_heterocyles": 0.0,
            "n_saturated_rings": 0.0,
        }
    )

    # Scenario #1
    props = dm.descriptors.compute_many_descriptors(mol)
    props = pd.Series(props)

    assert props.equals(true_values)

    # Scenario #2
    props = dm.descriptors.compute_many_descriptors(
        mol,
        properties_fn={"hello": lambda x: 88},
        add_properties=False,
    )
    assert props == {"hello": 88}

    # Scenario #3
    props = dm.descriptors.compute_many_descriptors(
        mol,
        properties_fn={"hello": lambda x: 88},
        add_properties=True,
    )
    props = pd.Series(props)

    true_values_2 = true_values.copy()
    true_values_2["hello"] = 88
    true_values_2 = true_values_2[props.index]

    assert true_values_2.equals(props)


def test_compute_many_descriptors_with_function_as_string():
    mol = dm.to_mol("CC(=O)OC1=CC=CC=C1C(=O)O")

    results = dm.descriptors.compute_many_descriptors(
        mol,
        properties_fn={"max_partial_charge": "MaxPartialCharge"},
        add_properties=False,
    )

    assert "max_partial_charge" in results.keys()
    assert pytest.approx(0.33900378687731025) == results["max_partial_charge"]


def test_batch_compute_many_descriptors():
    data = dm.data.freesolv()
    mols = data["smiles"].apply(dm.to_mol).tolist()

    props = dm.descriptors.batch_compute_many_descriptors(
        mols,
        batch_size=64,
        n_jobs=-1,
        progress=False,
    )

    assert set(props.columns.tolist()) == {
        "mw",
        "fsp3",
        "n_lipinski_hba",
        "n_lipinski_hbd",
        "n_rings",
        "n_hetero_atoms",
        "n_heavy_atoms",
        "n_rotatable_bonds",
        "n_radical_electrons",
        "tpsa",
        "qed",
        "clogp",
        "sas",
        "n_aliphatic_carbocycles",
        "n_aliphatic_heterocyles",
        "n_aliphatic_rings",
        "n_aromatic_carbocycles",
        "n_aromatic_heterocyles",
        "n_aromatic_rings",
        "n_saturated_carbocycles",
        "n_saturated_heterocyles",
        "n_saturated_rings",
    }
    assert props.shape == (642, 22)


def test_any_rdkit_descriptor():
    mol = dm.to_mol("CC(=O)OC1=CC=CC=C1C(=O)O")

    value = dm.descriptors.any_rdkit_descriptor("MaxPartialCharge")(mol)
    assert pytest.approx(value) == 0.33900378687731025

    value = dm.descriptors.any_rdkit_descriptor("CalcFractionCSP3")(mol)
    assert pytest.approx(value) == 0.1111111111111111

    with pytest.raises(ValueError):
        dm.descriptors.any_rdkit_descriptor("DOES NOT EXIST")


def test_n_aromatic_atoms():

    smiles = "Nc1cnn(-c2ccccc2)c(=O)c1Cl"
    mol = dm.to_mol(smiles)

    assert dm.descriptors.n_aromatic_atoms(mol) == 12
    assert dm.descriptors.n_aromatic_atoms_proportion(mol) == 0.8


def test_esol():

    smiles = "Nc1cnn(-c2ccccc2)c(=O)c1Cl"
    mol = dm.to_mol(smiles)

    assert np.allclose(dm.descriptors.esol(mol), -2.627091966265316)


def test_esol_from_data():
    data = dm.freesolv()
    data = data.iloc[:20]

    with pytest.raises(KeyError):
        dm.descriptors.esol_from_data(data)

    data["mol"] = data["smiles"].apply(dm.to_mol)
    data["clogp"] = data["mol"].apply(dm.descriptors.clogp)
    data["mw"] = data["mol"].apply(dm.descriptors.mw)
    data["n_rotatable_bonds"] = data["mol"].apply(dm.descriptors.n_rotatable_bonds)
    data["n_aromatic_atoms_proportion"] = data["mol"].apply(
        dm.descriptors.n_aromatic_atoms_proportion
    )

    # dataframe
    esol_values = dm.descriptors.esol_from_data(data)
    assert esol_values.dtype == float
    assert esol_values.shape == (20,)

    # series
    v = dm.descriptors.esol_from_data(data.iloc[0])
    v = float(v)
    assert type(v) == float

    # dict
    v = dm.descriptors.esol_from_data(data.iloc[0].to_dict())
    v = float(v)
    assert type(v) == float
