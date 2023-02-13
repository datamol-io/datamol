import pytest

import pandas as pd
import datamol as dm


def test_descriptors():
    smiles_list = ["CC(=O)OC1=CC=CC=C1C(=O)O", "CCN(CC)CCCC(C)NC1=C2C=CC(=CC2=NC=C1)Cl"]

    for smiles in smiles_list:
        mol = dm.to_mol(smiles)

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
        dm.descriptors.sas(mol)
        dm.descriptors.n_stereo_centers_unspecified(mol)
        dm.descriptors.n_spiro_atoms(mol)

        dm.descriptors.n_aliphatic_carbocycles(mol)
        dm.descriptors.n_aliphatic_heterocyles(mol)
        dm.descriptors.n_aliphatic_rings(mol)
        dm.descriptors.n_aromatic_carbocycles(mol)
        dm.descriptors.n_aromatic_heterocyles(mol)
        dm.descriptors.n_aromatic_rings(mol)
        dm.descriptors.n_saturated_carbocycles(mol)
        dm.descriptors.n_saturated_heterocyles(mol)
        dm.descriptors.n_saturated_rings(mol)


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
    data = data.iloc[:30]
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
    assert props.shape == (30, 22)


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


def test_formal_charge():
    mol = dm.to_mol("CC(=O)NC1CCC2=CC(=C(C(=C2C3=CC=C(C(=O)C=C13)OC)OC)OC)OC")
    assert dm.descriptors.formal_charge(mol) == 0

    mol = dm.to_mol("C(CC(=O)[O-])C(C(=O)[O-])[NH3+]")
    assert dm.descriptors.formal_charge(mol) == -1


def test_refractivity():
    mol = dm.to_mol("CN1C(=O)CN=C(C2=C1C=CC(=C2)Cl)C3=CC=CC=C3")

    value = dm.descriptors.refractivity(mol)
    assert pytest.approx(value, rel=2) == 81.10


def test_n_rigid_bonds():
    mol = dm.to_mol("CC(=O)NC1CCC2=CC(=C(C(=C2C3=CC=C(C(=O)C=C13)OC)OC)OC)OC")
    assert dm.descriptors.n_rigid_bonds(mol) == 20

    mol = dm.to_mol("CN1C(=O)CN=C(C2=C1C=CC(=C2)Cl)C3=CC=CC=C3")
    assert dm.descriptors.n_rigid_bonds(mol) == 19


def test_n_stereocenters():
    mol = dm.to_mol("CC(=O)NC1CCC2=CC(=C(C(=C2C3=CC=C(C(=O)C=C13)OC)OC)OC)OC")

    assert dm.descriptors.n_stereo_centers(mol) == 1

    mol = dm.to_mol("CN1C(=O)CN=C(C2=C1C=CC(=C2)Cl)C3=CC=CC=C3")
    assert dm.descriptors.n_stereo_centers(mol) == 0


def test_n_charged_atoms():
    mol = dm.to_mol("C(CC(=O)[O-])C(C(=O)[O-])[NH3+]")
    assert dm.descriptors.n_charged_atoms(mol) == 3
