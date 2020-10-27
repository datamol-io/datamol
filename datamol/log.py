import contextlib
import inspect

from rdkit import RDLogger


@contextlib.contextmanager
def without_rdkit_log(level=RDLogger.CRITICAL):
    """Context to disable rdkit log."""
    rdkit_logger = RDLogger.logger()
    rdkit_logger.setLevel(level)

    try:
        yield
    finally:
        rdkit_logger.setLevel(RDLogger.ERROR)


def disable_rdkit_log(level=RDLogger.CRITICAL):
    """Context to disable rdkit log."""
    rdkit_logger = RDLogger.logger()
    rdkit_logger.setLevel(level)
