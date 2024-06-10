from rdkit import RDLogger
from rdkit import rdBase
from functools import wraps


class without_rdkit_log:
    """Context manager to disable RDKit logs. By default all logs are disabled.

    Example:

    ```python
    import datamol as dm

    with dm.without_rdkit_log():
        mol = dm.to_mol("CCCCO")  # potential RDKit logs won't show
    ```
    """

    def __init__(
        self,
        mute_errors: bool = True,
        mute_warning: bool = True,
        mute_info: bool = True,
        mute_debug: bool = True,
        enable: bool = True,
    ):
        if enable is False:
            mute_errors = False
            mute_warning = False
            mute_info = False
            mute_debug = False

        # Get current log state
        self.previous_status = self._get_log_status()

        # Init the desired log state to apply during in the context
        self.desired_status = {}
        self.desired_status["rdApp.error"] = not mute_errors
        self.desired_status["rdApp.warning"] = not mute_warning
        self.desired_status["rdApp.debug"] = not mute_debug
        self.desired_status["rdApp.info"] = not mute_info

    def _get_log_status(self):
        """Get the current log status of RDKit logs."""
        log_status = rdBase.LogStatus()
        log_status = {st.split(":")[0]: st.split(":")[1] for st in log_status.split("\n")}
        log_status = {k: True if v == "enabled" else False for k, v in log_status.items()}
        return log_status

    def _apply_log_status(self, log_status):
        """Apply an RDKit log status."""
        for k, v in log_status.items():
            if v is True:
                rdBase.EnableLog(k)
            else:
                rdBase.DisableLog(k)

    def __enter__(self):
        self._apply_log_status(self.desired_status)

    def __exit__(self, *args, **kwargs):
        self._apply_log_status(self.previous_status)


def disable_rdkit_log():
    """Disable all rdkit logs."""
    for log_level in RDLogger._levels:
        rdBase.DisableLog(log_level)


def enable_rdkit_log():
    """Enable all rdkit logs."""
    for log_level in RDLogger._levels:
        rdBase.EnableLog(log_level)


def no_rdkit_log(
    func=None,
    *,
    mute_errors: bool = True,
    mute_warning: bool = True,
    mute_info: bool = True,
    mute_debug: bool = True,
    enable: bool = True,
):
    """Decorator to disable RDKit logs.

    This decorator can be used to suppress RDKit logs when executing a specific function.
    By default, all log levels (error, warning, info, and debug) are muted.

    Args:
        mute_errors : Whether to mute error logs (default is True).
        mute_warning : Whether to mute warning logs (default is True).
        mute_info : Whether to mute info logs (default is True).
        mute_debug : Whether to mute debug logs (default is True).
        enable: Whether to enable the log muting (default is True). If set to False, no logs will be muted.

    Example:
    ```python
    @no_rdkit_log()
    def example_function():
        # Your function code here
        pass

    example_function()  # RDKit logs won't show during this function's execution
    ```
    """

    if func is None:
        return lambda f: no_rdkit_log(
            f,
            mute_errors=mute_errors,
            mute_warning=mute_warning,
            mute_info=mute_info,
            mute_debug=mute_debug,
            enable=enable,
        )

    @wraps(func)
    def wrapper(*args, **kwargs):
        with without_rdkit_log(mute_errors, mute_warning, mute_info, mute_debug, enable):
            return func(*args, **kwargs)

    return wrapper
