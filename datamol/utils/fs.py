from typing import Union

import os

import fsspec


def copy_file(
    source: Union[str, os.PathLike],
    destination: Union[str, os.PathLike],
):
    """Copy a file to another path using fsspec.

    Args:
        source: Path to a file to copy from (remote or local).
        destination: Path to a file to copy to (remote or local).
    """
    with fsspec.open(destination, "wb") as f_dest:
        with fsspec.open(source, "rb") as f_source:
            f_dest.write(f_source.read())
