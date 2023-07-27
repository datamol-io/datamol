import platform
import pathlib
from loguru import logger

import pytest


DATA_DIR_PATH = pathlib.Path(__file__).parent.resolve() / "data"


@pytest.fixture
def current_platform():
    if platform.system() == "Linux":
        return "linux"
    elif platform.system() == "Darwin":
        return "osx"
    elif platform.system() == "Windows":
        return "win"
    else:
        return platform.system()


@pytest.fixture(autouse=True)
def skip_by_platform(request, current_platform):
    if request.node.get_closest_marker("skip_platform"):
        if request.node.get_closest_marker("skip_platform").args[0] == current_platform:
            pytest.skip(f"skipped on this platform: {current_platform}")


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "skip_platform(current_platform): skip test for a given platform from `['linux', 'osx', 'win']`",
    )


@pytest.fixture
def datadir(request):
    return DATA_DIR_PATH


# Mandatory for the below monkeypatch function.
from _pytest.logging import caplog as _caplog  # noqa: E402, F401


@pytest.fixture
def caplog(_caplog):  # noqa: F811
    """Monkeypatching the pytest caplog to work with loguru.

    See https://loguru.readthedocs.io/en/latest/resources/migration.html#making-things-work-with-pytest-and-caplog
    """
    import logging

    class PropogateHandler(logging.Handler):
        def emit(self, record):
            logging.getLogger(record.name).handle(record)

    handler_id = logger.add(PropogateHandler(), format="{message}")
    yield _caplog
    logger.remove(handler_id)
