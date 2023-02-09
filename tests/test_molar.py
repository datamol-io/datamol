import pytest

import datamol as dm
import numpy as np
import pandas as pd


MOLAR_TEST_VALUES = pd.DataFrame(
    [
        (1, 6, "uM"),
        (0.059, 7.229147988357856, "uM"),
        (0.024, 7.61978876, "uM"),
        (0.187, 6.72815839, "uM"),
        (0.00154, 8.8124793, "uM"),
        (128, 6.892790, "nM"),
        (0.000128, 6.892790, "mM"),
    ],
    columns=["xc50", "pxc50", "unit"],
)


def test_molar_to_log():
    # test scalar
    value, log_value, unit = MOLAR_TEST_VALUES.iloc[0].values
    assert dm.molar.molar_to_log(value, unit=unit) == log_value

    # test arrays
    for unit in ["uM", "mM", "nM"]:
        mask = MOLAR_TEST_VALUES["unit"] == unit
        values = MOLAR_TEST_VALUES[mask]["xc50"].tolist()
        log_values = MOLAR_TEST_VALUES[mask]["pxc50"].tolist()
        np.testing.assert_almost_equal(dm.molar.molar_to_log(values, unit=unit), log_values)

    # test wrong unit
    with pytest.raises(ValueError):
        dm.molar.molar_to_log(0.000128, unit="kcal/mol")


def test_log_to_molar():
    # test scalar
    value, log_value, unit = MOLAR_TEST_VALUES.iloc[0].values
    np.testing.assert_almost_equal(dm.molar.log_to_molar(log_value, unit=unit), value)

    # test arrays
    for unit in ["uM", "mM", "nM"]:
        mask = MOLAR_TEST_VALUES["unit"] == unit
        values = MOLAR_TEST_VALUES[mask]["xc50"].tolist()
        log_values = MOLAR_TEST_VALUES[mask]["pxc50"].tolist()
        np.testing.assert_almost_equal(
            dm.molar.log_to_molar(log_values, unit=unit), values, decimal=5
        )

    # test wrong unit
    with pytest.raises(ValueError):
        dm.molar.log_to_molar(7.214, unit="kcal/mol")


def test_log_to_molar_with_integer():
    dm.molar.log_to_molar(6, unit="uM")
