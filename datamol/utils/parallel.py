from typing import Iterable
from typing import Any
from typing import Callable
from typing import Optional

from concurrent import futures

from loguru import logger


def parallelized(
    fn: Callable,
    inputs_list: Iterable[Any],
    scheduler: str = "processes",
    n_jobs: Optional[int] = None,
    progress: bool = False,
    arg_type: str = "arg",
    _progress_auto: bool = True,
    _progress_leave: bool = False,
):
    """Run a function in parallel.

    Args:
        fn (Callable): The function to run in parallel.
        inputs_list (List[Any]): List of inputs to pass to `fn`.
        scheduler (str, optional): Choose between ["processes", "threads"]. Defaults
            to "processes".
        n_jobs (Optional[int], optional): Number of workers. If None, it will default
            to the number of processors on the machine. Defaults to None.
        progress (bool, optional): Display a progress bar. Defaults to False.
        arg_type (str, optional): One of ["arg", "args", "kwargs]:
            - "arg": the input is passed as an argument: `fn(arg)` (default).
            - "args": the input is passed as a list: `fn(*args)`.
            - "kwargs": the input is passed as a map: `fn(**kwargs)`.
        _progress_auto (bool, optional): Whether to use `tqdm.auto`. Defaults to True.
        _progress_leave (bool, optional): Leave progress bar after the execition. Defaults to False.

    Returns:
        The results of the execution as a list.
    """

    arg_type_values = ["arg", "args", "kwargs"]
    if not arg_type in arg_type_values:
        raise ValueError(f"`args_as` must be from {arg_type_values} instead of '{arg_type}'")

    if _progress_auto:
        from tqdm.auto import tqdm
    else:
        from tqdm import tqdm

    if scheduler == "processes":
        pool_executor_cls = futures.ProcessPoolExecutor
        executor_kwargs = dict(max_workers=n_jobs)
        as_completed_fn = futures.as_completed

    elif scheduler == "threads":
        pool_executor_cls = futures.ThreadPoolExecutor
        executor_kwargs = dict(max_workers=n_jobs)
        as_completed_fn = futures.as_completed

    else:
        raise ValueError(f"Wrong scheduler: {scheduler}")

    with pool_executor_cls(**executor_kwargs) as executor:

        if arg_type == "arg":
            futures_list = {executor.submit(fn, kwargs): i for i, kwargs in enumerate(inputs_list)}
        elif arg_type == "args":
            futures_list = {executor.submit(fn, *kwargs): i for i, kwargs in enumerate(inputs_list)}
        elif arg_type == "kwargs":
            futures_list = {
                executor.submit(fn, **kwargs): i for i, kwargs in enumerate(inputs_list)
            }
        else:
            raise ValueError("The following should never happen.")

        tqdm_args = {}
        try:
            tqdm_args["total"] = len(inputs_list)
        except Exception:
            tqdm_args["total"] = None
        tqdm_args["disable"] = not progress
        tqdm_args["leave"] = _progress_leave

        results = []
        for future in tqdm(as_completed_fn(futures_list), **tqdm_args):

            i = futures_list[future]
            try:
                result = future.result()
                results.append((i, result))
            except Exception as ex:
                logger.info("One worker failed.", exc_info=1)
                logger.info("Shutting down the execution...")
                executor.shutdown(True)
                raise ex

    # Reorder the results so they match the initial kwargs_list
    # NOTE(hadim): not sure about the performance here.
    results = sorted(results, key=lambda x: x[0])
    results = [result[1] for result in results]

    return results
