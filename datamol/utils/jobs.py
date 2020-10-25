from typing import Iterable
from typing import Callable
from typing import Optional
from typing import Any

from joblib import Parallel, delayed
from tqdm.auto import tqdm


class JobRunner:
    """
    JobRunner with sequential/parallel regimes. The multiprocessing backend use joblib which
    allows taking advantage of its features, while the progress bar use tqdm

    Args:
        n_jobs (int, optional): number of process (see joblib). Use 0 or None to force sequential.
        prefer (int, optional): Choose from ['processes', 'threads'] or None. Default to None.
            Soft hint to choose the default backend if no specific backend
            was selected with the parallel_backend context manager. The
            default process-based backend is 'loky' and the default
            thread-based backend is 'threading'. Ignored if the ``backend``
            parameter is specified.
        progress (bool, optional): whether to display progress bar
        job_kwargs (dict, optional): Any additional keyword argument supported by joblib.Parallel.

    Example:

        .. code-block:: python

            import datamol as dm
            runner = dm.JobRunner(n_jobs=4, progress=True, prefer="threads")
            results = runner(lambda x: x**2, [1, 2, 3, 4])

    """

    def __init__(self, n_jobs: int = 0, prefer: str = None, progress: bool = False, **job_kwargs):
        self.n_jobs = n_jobs
        self.prefer = prefer
        self.job_kwargs = job_kwargs
        self.job_kwargs.update(n_jobs=self.n_jobs, prefer=self.prefer)
        self.no_progress = not progress

    @property
    def is_sequential(self):
        """Check whether the job is sequential or parallel"""
        return (self.n_jobs is None) or (self.n_jobs in [0, 1])

    @staticmethod
    def wrap_fn(fn: Callable, arg_type: Optional[str] = None, **fn_kwargs):
        """Small wrapper around a callable to properly format it's argument"""
        # EN probably use something like (moms.utils.commons.is_callable) ?
        def _run(args: Any):
            if arg_type == "kwargs":
                fn_kwargs.update(**args)
                return fn(**fn_kwargs)
            elif arg_type == "args":
                return fn(*args, **fn_kwargs)
            return fn(args, **fn_kwargs)

        return _run

    @staticmethod
    def get_iterator_length(data):
        """Attempt to get the length of an iterator"""
        total_length = None
        try:
            total_length = len(data)
        except TypeError:
            # most likely a generator, ignore
            pass
        return total_length

    def sequential(
        self,
        callable_fn: Callable,
        data: Iterable[Any],
        arg_type: Optional[str] = None,
        **fn_kwargs,
    ):
        r"""
        Run job in sequential version

        Args:
            callable_fn (callable): function to call
            data (iterable): input data
            arg_type (str, optional): function argument type ('arg'/None or 'args' or 'kwargs')
            fn_kwargs (dict, optional): optional keyword argument to pass to the callable funciton
        """
        total_length = JobRunner.get_iterator_length(data)
        res = [
            JobRunner.wrap_fn(callable_fn, arg_type, **fn_kwargs)(dt)
            for dt in tqdm(data, total=total_length, disable=self.no_progress)
        ]
        return res

    def parallel(
        self,
        callable_fn: Callable,
        data: Iterable[Any],
        arg_type: Optional[str] = None,
        **fn_kwargs,
    ):
        r"""
        Run job in parallel

        Args:
            callable_fn (callable): function to call
            data (iterable): input data
            arg_type (str, optional): function argument type ('arg'/None or 'args' or 'kwargs')
            fn_kwargs (dict, optional): optional keyword argument to pass to the callable funciton
        """
        runner = JobRunner._parallel_helper(**self.job_kwargs)
        total_length = JobRunner.get_iterator_length(data)
        results = runner(total=total_length, disable=self.no_progress)(
            delayed(JobRunner.wrap_fn(callable_fn, arg_type, **fn_kwargs))(dt) for dt in data
        )
        return results

    def __call__(self, *args, **kwargs):
        """
        Run job using the n_jobs attribute to determine regime
        """
        if self.is_sequential:
            return self.sequential(*args, **kwargs)
        return self.parallel(*args, **kwargs)

    @staticmethod
    def _parallel_helper(**joblib_args):
        r"""
        Parallel helper function for joblib with tqdm support
        """

        def run(**tq_args):
            def tmp(op_iter):
                return Parallel(**joblib_args)(tqdm(op_iter, **tq_args))

            return tmp

        return run
