import pytest

import pandas as pd
import datamol as dm


def test_descriptors():

    mol = dm.to_mol("CCN(CC)CCCC(C)NC1=C2C=CC(=CC2=NC=C1)Cl")
    dm.descriptors.mw(mol)
    dm.descriptors.fsp3(mol)
    dm.descriptors.n_hba(mol)
    dm.descriptors.n_hbd(mol)
    dm.descriptors.n_rings(mol)
    dm.descriptors.n_hetero_atoms(mol)
    dm.descriptors.n_heavy_atoms(mol)
    dm.descriptors.n_rotatable_bonds(mol)
    dm.descriptors.n_aliphatic_rings(mol)
    dm.descriptors.n_aromatic_rings(mol)
    dm.descriptors.n_saturated_rings(mol)
    dm.descriptors.n_H_acceptors(mol)
    dm.descriptors.n_H_donors(mol)
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
            "n_hba": 3,
            "n_hbd": 1,
            "n_rings": 2,
            "n_hetero_atoms": 4,
            "n_heavy_atoms": 22,
            "n_rotatable_bonds": 8,
            "n_aliphatic_rings": 0,
            "n_aromatic_rings": 2,
            "n_saturated_rings": 0,
            "n_H_acceptors": 3,
            "n_H_donors": 1,
            "n_radical_electrons": 0,
            "tpsa": 28.16,
            "qed": 0.7564117572128701,
            "clogp": 4.810600000000004,
            "sas": 2.670786229594949,
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


def test_batch_compute_many_descriptors():
    data = dm.data.freesolv()
    mols = data["smiles"].apply(dm.to_mol).tolist()

    props = dm.descriptors.batch_compute_many_descriptors(
        mols,
        batch_size=64,
        n_jobs=-1,
        progress=False,
    )

    assert props.columns.tolist() == [
        "mw",
        "fsp3",
        "n_hba",
        "n_hbd",
        "n_rings",
        "n_hetero_atoms",
        "n_heavy_atoms",
        "n_rotatable_bonds",
        "n_aliphatic_rings",
        "n_aromatic_rings",
        "n_saturated_rings",
        "n_H_acceptors",
        "n_H_donors",
        "n_radical_electrons",
        "tpsa",
        "qed",
        "clogp",
        "sas",
    ]
    assert props.shape == (642, 18)


def test_any_descriptor():
    mol = dm.to_mol("CC(=O)OC1=CC=CC=C1C(=O)O")

    value = dm.descriptors.any_descriptor("MaxPartialCharge")(mol)
    assert pytest.approx(value) == 0.33900378687731025

    value = dm.descriptors.any_descriptor("CalcFractionCSP3")(mol)
    assert pytest.approx(value) == 0.1111111111111111

    with pytest.raises(ValueError):
        dm.descriptors.any_descriptor("DOES NOT EXIST")
