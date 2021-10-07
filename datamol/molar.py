"""A set of utility functions to convert between various units and formats used in drug discovery.
"""

from typing import Union
from typing import Iterable

import numpy as np


_MOLAR_SCALES = {"M": 1, "mM": 1e-3, "uM": 1e-6, "nM": 1e-9, "pM": 1e-12, "fM": 1e-15}


def molar_to_log(
    values: Union[float, Iterable[float], np.ndarray],
    unit: str,
) -> Union[float, Iterable[float], np.ndarray]:
    """Convert a molar concentration (XC50 for example) to its log scaled value (pXC50).

    Args:
        values: A molar concentration (can be a scalar, a list or an array).
        unit: The unit of the input concentration. Choose from:
            `{'M', 'fM', 'mM', 'nM', 'pM', 'uM'}`.
    """

    if unit not in _MOLAR_SCALES:
        raise ValueError(
            f"The unit '{unit}' is not supported. Choose from {set(_MOLAR_SCALES.keys())}."
        )

    return -1 * np.log10(np.array(values) * _MOLAR_SCALES[unit])


def log_to_molar(
    values: Union[float, Iterable[float], np.ndarray],
    unit: str,
) -> Union[float, Iterable[float], np.ndarray]:
    """Convert a log-scaled molar concentration (pXC50 for example) to its unscaled value (XC50).

    Args:
        values: A log-scaled molar concentration (can be a scalar, a list or an array).
        unit: The unit of the input concentration. Choose from:
            `{'M', 'fM', 'mM', 'nM', 'pM', 'uM'}`.
    """

    if unit not in _MOLAR_SCALES:
        raise ValueError(
            f"The unit '{unit}' is not supported. Choose from {set(_MOLAR_SCALES.keys())}."
        )

    return 10 ** (-1 * np.array(values, dtype="float")) / _MOLAR_SCALES[unit]
