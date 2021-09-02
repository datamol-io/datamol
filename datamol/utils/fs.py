"""The `fs` module makes it easier to work with all type of path (the ones supported by `fsspec`).
"""

from typing import Union
from typing import Optional
from typing import List

import os
import io
import hashlib
import pathlib

import fsspec


def _import_tqdm():
    try:
        from tqdm.auto import tqdm

        return tqdm
    except ImportError:
        return None


def _import_appdirs():
    try:
        import appdirs

        return appdirs
    except ImportError:
        return None


def get_cache_dir(app_name: str, suffix: str = None, create: bool = True):
    """Get a local cache directory for a given application name.

    Args:
        app_name: The name of the application.
        suffix: A subdirectory appended to the cache dir.
        create: Whether to create the directory and its parents if it does not
            already exist.
    """

    appdirs = _import_appdirs()

    if appdirs is None:
        raise ImportError(
            "To use `dm.utils.fs.get_cache_dir()`, you must have `appdirs` "
            "installed: `conda install appdirs`."
        )
    cache_dir = pathlib.Path(appdirs.user_cache_dir(appname=app_name))

    if suffix is not None:
        cache_dir /= suffix

    if create:
        cache_dir.mkdir(exist_ok=True, parents=True)

    return cache_dir


def get_mapper(path: Union[str, os.PathLike]):
    """Get the fsspec mapper.

    Args:
        path: a path supported by `fsspec` such as local, s3, gcs, etc.
    """
    return fsspec.get_mapper(str(path))


def get_basename(path: Union[str, os.PathLike]):
    """Get the basename of a file or a folder.

    Args:
        path: a path supported by `fsspec` such as local, s3, gcs, etc.
    """
    path = str(path)
    mapper = get_mapper(path)
    clean_path = path.rstrip(mapper.fs.sep)
    return str(clean_path).split(mapper.fs.sep)[-1]


def get_extension(path: Union[str, os.PathLike]):
    """Get the extension of a file.

    Args:
        path: a path supported by `fsspec` such as local, s3, gcs, etc.
    """
    basename = get_basename(path)
    return basename.split(".")[-1]


def exists(path: Union[str, os.PathLike, fsspec.core.OpenFile, io.IOBase]):
    """Check whether a file or a directory exists.

    Important: File-like object always exists.

    Args:
        path: a path supported by `fsspec` such as local, s3, gcs, etc.
    """
    return is_file(path) or is_dir(path)


def is_file(path: Union[str, os.PathLike, fsspec.core.OpenFile, io.IOBase]):
    """Check whether a file exists.

    Important: File-like object always exists.

    Args:
        path: a path supported by `fsspec` such as local, s3, gcs, etc.
    """
    if isinstance(path, fsspec.core.OpenFile):
        return path.fs.isfile(path.path)

    elif isinstance(path, (str, pathlib.Path)):
        mapper = get_mapper(str(path))
        return mapper.fs.isfile(path)

    else:
        return True


def is_dir(path: Union[str, os.PathLike, fsspec.core.OpenFile, io.IOBase]):
    """Check whether a file exists.

    Important: File-like object always exists.

    Args:
        path: a path supported by `fsspec` such as local, s3, gcs, etc.
    """
    if isinstance(path, fsspec.core.OpenFile):
        return path.fs.isdir(path.path)

    elif isinstance(path, (str, pathlib.Path)):
        mapper = get_mapper(str(path))
        return mapper.fs.isdir(path)

    else:
        return False


def get_protocol(path: Union[str, os.PathLike]):
    """Return the name of the path protocol.

    Args:
        path: a path supported by `fsspec` such as local, s3, gcs, etc.
    """

    mapper = get_mapper(path)
    protocol = mapper.fs.protocol

    if "s3" in protocol:
        return "s3"
    elif "gs" in protocol:
        return "gs"
    elif isinstance(protocol, (tuple, list)):
        return protocol[0]
    return protocol


def is_local_path(path: Union[str, os.PathLike]):
    """Check whether a path is local."""
    return get_protocol(str(path)) == "file"


def join(*paths):
    """Join paths together. The first element determine the
    filesystem to use (and so the separator.

    Args:
        paths: a list of paths supported by `fsspec` such as local, s3, gcs, etc.
    """
    paths = [str(path).rstrip("/") for path in paths]
    source_path = paths[0]
    fs = get_mapper(source_path).fs
    full_path = fs.sep.join(paths)
    return full_path


def get_size(file: Union[str, os.PathLike, io.IOBase, fsspec.core.OpenFile]) -> Optional[int]:
    """Get the size of a file given its path. Return None if the
    size can't be retrieved.
    """

    if isinstance(file, io.IOBase) and hasattr(file, "name"):
        fs_local = fsspec.filesystem("file")
        file_size = fs_local.size(getattr(file, "name"))

    elif isinstance(file, (str, pathlib.Path)):
        fs = get_mapper(str(file)).fs
        file_size = fs.size(str(file))

    elif isinstance(file, fsspec.core.OpenFile):
        file_size = file.fs.size(file.path)

    else:
        file_size = None

    return file_size


def copy_file(
    source: Union[str, pathlib.Path, io.IOBase, fsspec.core.OpenFile],
    destination: Union[str, pathlib.Path, io.IOBase, fsspec.core.OpenFile],
    chunk_size: int = None,
    force: bool = False,
    progress: bool = False,
    leave_progress: bool = True,
):
    """Copy one file to another location across different filesystem (local, S3, GCS, etc).

    Args:
        source: path or file-like object to copy from.
        destination: path or file-like object to copy to.
        chunk_size: the chunk size to use. If progress is enabled the chunk
            size is `None`, it is set to 2048.
        force: whether to overwrite the destination file it it exists.
        progress: whether to display a progress bar.
        leave_progress: whether to hide the progress bar once the copy is done.
    """

    if progress and chunk_size is None:
        chunk_size = 2048

    if isinstance(source, (str, pathlib.Path)):
        source_file = fsspec.open(str(source), "rb")
    else:
        source_file = source

    if isinstance(destination, (str, pathlib.Path)):

        # adapt the file mode of the destination depending on the source file.
        destination_mode = "wb"
        if hasattr(source_file, "mode"):
            destination_mode = "wb" if "b" in getattr(source_file, "mode") else "w"
        elif isinstance(source_file, io.BytesIO):
            destination_mode = "wb"
        elif isinstance(source_file, io.StringIO):
            destination_mode = "w"

        destination_file = fsspec.open(str(destination), destination_mode)
    else:
        destination_file = destination

    if not is_file(source_file):
        raise ValueError(f"The file being copied does not exist: {source}")

    if not force and is_file(destination_file):
        raise ValueError(f"The destination file to copy already exists: {destination}")

    with source_file as source_stream:
        with destination_file as destination_stream:

            if chunk_size is None:
                # copy without chunks
                destination_stream.write(source_stream.read())

            else:
                # copy with chunks

                # determine the size of the source file
                source_size = None
                if progress:
                    source_size = get_size(source)

                pbar = None
                if progress:
                    tqdm = _import_tqdm()

                    if tqdm is None:
                        raise ImportError(
                            "If the progress bar is enabled, you must have `tqdm` "
                            "installed: `conda install tqdm`."
                        )
                    else:
                        # init progress bar
                        pbar = tqdm(
                            total=source_size,
                            leave=leave_progress,
                            disable=not progress,
                            unit="B",
                            unit_divisor=1024,
                            unit_scale=True,
                        )

                # start the loop
                while True:
                    data = source_stream.read(chunk_size)
                    if not data:
                        break
                    destination_stream.write(data)

                    if pbar is not None:
                        pbar.update(chunk_size)

                if pbar is not None:
                    pbar.close()


def md5(filepath: Union[str, os.PathLike]):
    """Return the md5 hash of a file.

    NOTE(hadim): Use fsspec caching here maybe.

    Args:
        filepath: The path to the file to compute the MD5 hash on.
    """
    with fsspec.open(filepath) as f:
        file_hash = hashlib.md5()
        file_hash.update(f.read())
        file_hash = file_hash.hexdigest()
    return file_hash


def glob(path: str, **kwargs) -> List[str]:
    """Find files by glob-matching.

    Args:
        path: A glob-style path.
    """
    # Get the list of paths
    fs = get_mapper(path).fs
    data_paths = fs.glob(path, **kwargs)
    protocol = get_protocol(path)
    # Append path prefix if needed
    if protocol not in ["file", "https", "http"]:
        data_paths = [f"{protocol}://{d}" for d in data_paths]

    return data_paths
