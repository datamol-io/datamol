from typing import Optional

import time
import signal

from loguru import logger


duration_intervals = (
    ("weeks", 604800),  # 60 * 60 * 24 * 7
    ("days", 86400),  # 60 * 60 * 24
    ("h", 3600),  # 60 * 60
    ("min", 60),
    ("s", 1),
    ("ms", 1e-3),
    ("us", 1e-6),
)


def human_duration(seconds: float, granularity: int = 1):

    # NOTE(hadim): far from being perfect.

    result = []
    duration: float = seconds
    for name, count in duration_intervals:
        value = duration // count
        if value:
            duration -= value * count
            result.append(f"{value:.0f}{name}")
    return ", ".join(result[:granularity])


class watch_duration:
    """A Python decorator to measure execution time with logging capability.

    Args:
        log: Whether to log the measured duration.
        log_human_duration: Whether to log duration in a human way
            depending on the amount.

    Example:

    ```python
    def fn(n):
        for i in range(n):
            print(i)
            time.sleep(0.2)

    with dm.utils.perf.watch_duration(log=True) as w:
        fn(5)

    print(w.duration)
    ```
    """

    def __init__(self, log: bool = True, log_human_duration: bool = True):
        self.log = log
        self.log_human_duration = log_human_duration

        self.start = None
        self.end = None
        self.duration = None
        self.duration_minutes = None

    def __enter__(self):
        self.start = time.time()
        return self

    def __exit__(self, *_):

        assert self.start is not None

        self.end = time.time()
        self.duration = self.end - self.start
        self.duration_minutes = self.duration / 60

        if self.log:
            if self.log_human_duration:
                logger.info(f"Duration {human_duration(self.duration)}.")
            else:
                logger.info(f"Duration {self.duration_minutes:.2f} minutes")
