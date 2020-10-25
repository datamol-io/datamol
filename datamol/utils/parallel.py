from concurrent import futures

from typing import Iterable
from typing import Any
from typing import Callable
from typing import Optional

from loguru import logger

from .jobs import JobRunner


def parallelized(
    fn: Callable,
    inputs_list: Iterable[Any],
    scheduler: str = "processes",
    n_jobs: Optional[int] = None,
    progress: bool = False,
    arg_type: str = "arg",
):
    """Run a function in parallel.

    Args:
        fn (Callable): The function to run in parallel.
        inputs_list (List[Any]): List of inputs to pass to `fn`.
        scheduler (str, optional): Choose between ["processes", "threads"]. Defaults
            to None which uses the default joblib "loky" scheduler.
        n_jobs (Optional[int], optional): Number of workers. If None, it will default
            to the number of processors on the machine. Defaults to None.
        progress (bool, optional): Display a progress bar. Defaults to False.
        arg_type (str, optional): One of ["arg", "args", "kwargs]:
            - "arg": the input is passed as an argument: `fn(arg)` (default).
            - "args": the input is passed as a list: `fn(*args)`.
            - "kwargs": the input is passed as a map: `fn(**kwargs)`.
    Returns:
        The results of the execution as a list.
    """

    runner = JobRunner(n_jobs=n_jobs, progress=progress, prefer=scheduler)
    results = runner(fn, inputs_list, arg_type=arg_type)
    return results
