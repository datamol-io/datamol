import math
import numbers
import operator
import unittest

import numpy as np
import pandas as pd

from functools import reduce

import datamol as dm


def random_fn(*args, op="mul", **kwargs):
    """Perform random functions on a list"""
    all_values = [x for x in args if isinstance(x, numbers.Number)]
    all_values += [x for x in kwargs.values() if isinstance(x, numbers.Number)]
    op_fn = getattr(operator, op, None)
    if op_fn is None:
        op_fn = getattr(math, op)
        return op_fn(all_values[0])
    return reduce(op_fn, all_values)


class TestJobs(unittest.TestCase):
    def test_sequential(self):
        jobrunner = dm.JobRunner(n_jobs=None, progress=False)
        # practically do nothing (add a single value with nothing)
        o1 = jobrunner(random_fn, [9, 25, 1024], op="add")
        self.assertEqual(o1, [9, 25, 1024])

        # take the sqrt
        o2 = jobrunner(random_fn, [9, 25, 1024], op="sqrt")
        self.assertEqual(o2, [3, 5, 32])

        # multiply all inputs
        o3 = jobrunner(random_fn, [(1, 2, 3), (4, 5, 6), (3, 4, 0)], arg_type="args", op="mul")
        self.assertEqual(o3, [6, 4 * 5 * 6, 0])

        # do the same thing but with kwargs
        o4 = jobrunner(
            random_fn,
            iter([dict(a=1, b=2, c=3), dict(a=4, b=5, c=6), dict(a=3, b=4, c=0)]),
            arg_type="kwargs",
            op="mul",
        )
        self.assertEqual(o4, [6, 4 * 5 * 6, 0])

        o5 = jobrunner(random_fn, np.asarray([9, 25, 1024]), op="add")
        self.assertEqual(o5, [9, 25, 1024])

    def test_parallel(self):
        jobrunner1 = dm.JobRunner(n_jobs=4, progress=True)  # use loky backend
        o1 = jobrunner1(random_fn, [9, 25, 1024], op="add")
        self.assertEqual(o1, [9, 25, 1024])

        o5 = jobrunner1(random_fn, np.asarray([9, 25, 1024]), op="add")
        self.assertEqual(o5, [9, 25, 1024])

        o3 = jobrunner1(random_fn, [(1, 2, 3), (4, 5, 6), (3, 4, 0)], arg_type="args", op="mul")
        self.assertEqual(o3, [6, 4 * 5 * 6, 0])

        # use threads instead, no progress
        jobrunner2 = dm.JobRunner(n_jobs=2, progress=False, prefer="threads")
        o2 = jobrunner2(random_fn, [9, 25, 1024], op="sqrt")
        self.assertEqual(o2, [3, 5, 32])

        o4 = jobrunner2(
            random_fn,
            iter([dict(a=1, b=2, c=3), dict(a=4, b=5, c=6), dict(a=3, b=4, c=0)]),
            arg_type="kwargs",
            op="mul",
        )
        self.assertEqual(o4, [6, 4 * 5 * 6, 0])

    def test_seq_vs_parallel(self):
        # test parallel vs sequential
        jobrunner = dm.JobRunner(n_jobs=4, progress=False)  # use loky backend
        o_seq = jobrunner.sequential(
            random_fn, [(1, 2, 3), (4, 5, 6), (3, 4, 0)], arg_type="args", op="mul"
        )
        o_par = jobrunner.parallel(
            random_fn, [(1, 2, 3), (4, 5, 6), (3, 4, 0)], arg_type="args", op="mul"
        )
        self.assertEqual(o_seq, o_par)

    def test_parallelized(self):
        def fn(x):
            return x**2

        results = dm.parallelized(
            fn,
            [{"x": i} for i in range(10)],
            scheduler="processes",
            n_jobs=None,
            arg_type="kwargs",
            progress=True,
        )
        assert results == [0, 1, 4, 9, 16, 25, 36, 49, 64, 81]

        results = dm.parallelized(
            fn,
            [[i] for i in range(10)],
            scheduler="processes",
            n_jobs=None,
            arg_type="args",
            progress=True,
        )
        assert results == [0, 1, 4, 9, 16, 25, 36, 49, 64, 81]

        results = dm.parallelized(
            fn,
            range(10),
            scheduler="processes",
            n_jobs=None,
            progress=False,
        )
        assert results == [0, 1, 4, 9, 16, 25, 36, 49, 64, 81]

    def test_job_kwargs(self):
        def fn(x):
            return x**2

        results = dm.parallelized(
            fn,
            [{"x": i} for i in range(10)],
            scheduler="processes",
            n_jobs=None,
            arg_type="kwargs",
            progress=True,
            verbose=100,
        )
        assert results == [0, 1, 4, 9, 16, 25, 36, 49, 64, 81]

    def test_tqdm_kwargs(self):
        def fn(x):
            return x**2

        results = dm.parallelized(
            fn,
            [{"x": i} for i in range(10)],
            scheduler="processes",
            n_jobs=None,
            arg_type="kwargs",
            progress=True,
            tqdm_kwargs=dict(desc="My progress bar"),
        )
        assert results == [0, 1, 4, 9, 16, 25, 36, 49, 64, 81]

    def test_with_batch_size(self):
        def _fn(n):
            return n * 3

        def _fn_return_none(n):
            return None

        results = dm.utils.parallelized(
            _fn,
            range(997),
            n_jobs=-1,
            progress=True,
            batch_size=10,
        )
        assert len(results) == 997

        results = dm.utils.parallelized(
            _fn_return_none,
            range(997),
            n_jobs=-1,
            progress=True,
            batch_size=10,
        )
        assert len(results) == 997

    def test_with_total(self):
        def _fn_process_fn(_, row):
            datum = {}
            datum["smiles"] = row["smiles"]
            return pd.Series(datum)

        data = dm.freesolv()
        data = data.iloc[:50]

        # parallel mode

        ## check the `total` arg is ok
        dm.parallelized(
            _fn_process_fn,
            data.iterrows(),
            n_jobs=-1,
            progress=True,
            arg_type="args",
            total=50,
        )

        ## check collision between guessed total and provided one
        dm.parallelized(
            _fn_process_fn,
            list(data.iterrows()),
            n_jobs=-1,
            progress=True,
            arg_type="args",
            total=50,
        )

        # sequential mode

        ## check the `total` arg is ok
        dm.parallelized(
            _fn_process_fn,
            data.iterrows(),
            n_jobs=1,
            progress=True,
            arg_type="args",
            total=50,
        )

        ## check collision between guessed total and provided one
        dm.parallelized(
            _fn_process_fn,
            list(data.iterrows()),
            n_jobs=1,
            progress=True,
            arg_type="args",
            total=50,
        )


def test_parallelized_with_batches():
    data = dm.freesolv()
    data = data.iloc[:10]

    def _fn1(smiles):
        return len(smiles)

    results1 = dm.parallelized(
        _fn1,
        data["smiles"],
        progress=False,
        n_jobs=-1,
    )

    def _fn2(smiles_list):
        return [len(s) for s in smiles_list]

    results2 = dm.parallelized_with_batches(
        _fn2,
        data["smiles"],
        batch_size=2,
        progress=False,
        n_jobs=-1,
    )

    assert results1 == results2
