from rdkit import RDLogger
from rdkit import rdBase


class without_rdkit_log:
    """Context manager to disable RDKit logs. By default all logs are disabled."""

    def __init__(
        self,
        mute_errors: bool = True,
        mute_warning: bool = True,
        mute_info: bool = True,
        mute_debug: bool = True,
    ):
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
