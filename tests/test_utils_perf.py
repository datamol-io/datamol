import pytest

import time

import datamol as dm


def test_watch_duration():
    def fn(n):
        for i in range(n):
            print(i)

    with dm.utils.perf.watch_duration(log=True) as w:
        fn(5)

    assert isinstance(w.duration, float)


@pytest.mark.skip_platform("win")
def test_timeout():
    with dm.utils.perf.Timeout(seconds=2):
        value = None
        try:
            for i in range(5):
                value = i
                time.sleep(1.5)
        except TimeoutError:
            pass

    assert value == 1


@pytest.mark.skip_platform("win")
def test_timeout_none():
    with dm.utils.perf.Timeout(seconds=None):
        value = None
        try:
            for i in range(5):
                value = i
                time.sleep(0.05)
        except TimeoutError:
            pass

    assert value == 4
