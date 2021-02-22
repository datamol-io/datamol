from typing import Callable
from typing import List
from typing import Union

import platform
from functools import wraps


def disable_on_os(os_names: Union[str, List[str]]):
    """A decorator to disable a function raising an error if the OS detected is not supported.

    Args:
        os_names: OS names to disable this function. Valid OS names are: `["linux", "osx", "win"]`.
    """

    if isinstance(os_names, str):
        os_names = [os_names]

    valid_os_names = []
    for os_name in os_names:
        if os_name == "linux":
            valid_os_names.append("Linux")
        elif os_name == "win":
            valid_os_names.append("Windows")
        elif os_name == "osx":
            valid_os_names.append("Darwin")
        else:
            valid_os_names.append(os_name)

    def real_decorator(function: Callable):
        @wraps(function)
        def wrapper(*args, **kwargs):

            if platform.system() not in valid_os_names:
                retval = function(*args, **kwargs)
                return retval
            else:
                raise NotImplementedError(
                    f"The function {function.__name__} is not supported"
                    f" for the platform '{platform.system()}'."
                )

        return wrapper

    return real_decorator
