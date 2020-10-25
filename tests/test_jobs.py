import math
import numbers
import operator
import unittest
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

    def test_parallel(self):
        jobrunner1 = dm.JobRunner(n_jobs=4, progress=True)  # use loky backend
        o1 = jobrunner1(random_fn, [9, 25, 1024], op="add")
        self.assertEqual(o1, [9, 25, 1024])

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
        o_seq = jobrunner(random_fn, [(1, 2, 3), (4, 5, 6), (3, 4, 0)], arg_type="args", op="mul")
        o_par = jobrunner(random_fn, [(1, 2, 3), (4, 5, 6), (3, 4, 0)], arg_type="args", op="mul")
        self.assertEqual(o_seq, o_par)
