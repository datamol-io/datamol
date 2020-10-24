import unittest

import datamol as dm


def fn(x):
    return x ** 2


class TestParallel(unittest.TestCase):
    def test_parallelized(self):

        results = dm.parallelized(
            fn,
            [{"x": i} for i in range(10)],
            scheduler="processes",
            max_workers=None,
            arg_type="kwargs",
            progress=True,
        )
        assert results == [0, 1, 4, 9, 16, 25, 36, 49, 64, 81]

        results = dm.parallelized(
            fn,
            [[i] for i in range(10)],
            scheduler="processes",
            max_workers=None,
            arg_type="args",
            progress=True,
        )
        assert results == [0, 1, 4, 9, 16, 25, 36, 49, 64, 81]

        results = dm.parallelized(
            fn,
            [i for i in range(10)],
            scheduler="processes",
            max_workers=None,
            progress=False,
        )
        assert results == [0, 1, 4, 9, 16, 25, 36, 49, 64, 81]
