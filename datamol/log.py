import contextlib

from rdkit import RDLogger


@contextlib.contextmanager
def disable_rdkit_log(level=RDLogger.CRITICAL):
    """Disable rdkit log."""
    rdkit_logger = RDLogger.logger()
    rdkit_logger.setLevel(level)
    try:
        yield
    finally:
        rdkit_logger.setLevel(RDLogger.ERROR)
