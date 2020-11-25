from typing import Union
from typing import Optional
from typing import List
from typing import Sequence

import pathlib

import fsspec


def copy_files(
    source: Union[str, pathlib.Path],
    destination: Union[str, pathlib.Path],
):
    """Copy a file to another path using fsspec.

    Args:
        source: Path to the file. Path to the file. Can be a local path, HTTP,
            HTTPS, S3, GS, etc. See https://filesystem-spec.readthedocs.io/en/latest/
        destination: Path to the file. Path to the file. Can be a local path, HTTP,
            HTTPS, S3, GS, etc. See https://filesystem-spec.readthedocs.io/en/latest/
    """
    with fsspec.open(destination, "wb") as f_dest:
        with fsspec.open(source, "rb") as f_source:
            f_dest.write(f_source.read())
