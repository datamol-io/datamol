import platform

import pytest


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
        "skip_platform(current_platform): skip test for the given search engine",
    )
