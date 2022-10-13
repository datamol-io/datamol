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
import fsspec.utils

from .jobs import parallelized


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


def get_cache_dir(app_name: str, suffix: Optional[str] = None, create: bool = True):
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

    Args:
        path: a path supported by `fsspec` such as local, s3, gcs, etc.
    """
    if isinstance(path, fsspec.core.OpenFile):
        return path.fs.isfile(path.path)

    elif isinstance(path, (str, os.PathLike)):
        mapper = get_mapper(str(path))
        return mapper.fs.isfile(str(path))

    else:
        return False


def is_dir(path: Union[str, os.PathLike, fsspec.core.OpenFile, io.IOBase]):
    """Check whether a file exists.

    Args:
        path: a path supported by `fsspec` such as local, s3, gcs, etc.
    """
    if isinstance(path, fsspec.core.OpenFile):
        return path.fs.isdir(path.path)

    elif isinstance(path, (str, os.PathLike)):
        mapper = get_mapper(str(path))
        return mapper.fs.isdir(str(path))

    else:
        return False


def get_protocol(path: Union[str, os.PathLike], fs: Optional[fsspec.AbstractFileSystem] = None):
    """Return the name of the path protocol.

    Args:
        path: a path supported by `fsspec` such as local, s3, gcs, etc.
    """

    if fs is None:
        fs = get_mapper(path).fs

    protocol = fs.protocol  # type: ignore

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


def join(*paths: str):
    """Join paths together. The first element determine the
    filesystem to use (and so the separator.

    Args:
        *paths: a list of paths supported by `fsspec` such as local, s3, gcs, etc.
    """
    _paths = [str(path).rstrip("/") for path in paths]
    source_path = _paths[0]
    fs = get_mapper(source_path).fs
    full_path = fs.sep.join(_paths)
    return full_path


def get_size(file: Union[str, os.PathLike, io.IOBase, fsspec.core.OpenFile]) -> Optional[int]:
    """Get the size of a file given its path. Return None if the
    size can't be retrieved.
    """

    if isinstance(file, io.IOBase) and hasattr(file, "name"):
        fs_local = fsspec.filesystem("file")
        file_size = fs_local.size(getattr(file, "name"))

    elif isinstance(file, (str, os.PathLike)):
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
    chunk_size: Optional[int] = None,
    force: bool = False,
    progress: bool = False,
    leave_progress: bool = True,
):
    """Copy one file to another location across different filesystem (local, S3, GCS, etc).

    Args:
        source: path or file-like object to copy from.
        destination: path or file-like object to copy to.
        chunk_size: the chunk size to use. If progress is enabled the chunk
            size is `None`, it is set to 1MB (1024 * 1024).
        force: whether to overwrite the destination file if it exists.
        progress: whether to display a progress bar.
        leave_progress: whether to hide the progress bar once the copy is done.
    """

    if progress and chunk_size is None:
        chunk_size = 1024 * 1024

    if isinstance(source, (str, os.PathLike)):
        source_file = fsspec.open(str(source), "rb")
    else:
        source_file = source

    if isinstance(destination, (str, os.PathLike)):

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

    if not is_file(source_file):  # type: ignore
        raise ValueError(f"The file being copied does not exist or is not a file: {source}")

    if not force and is_file(destination_file):  # type: ignore
        raise ValueError(f"The destination file to copy already exists: {destination}")

    with source_file as source_stream:
        with destination_file as destination_stream:

            if chunk_size is None:
                # copy without chunks
                destination_stream.write(source_stream.read())  # type: ignore

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
                    data = source_stream.read(chunk_size)  # type: ignore
                    if not data:
                        break
                    destination_stream.write(data)  # type: ignore

                    if pbar is not None:
                        pbar.update(chunk_size)

                if pbar is not None:
                    pbar.close()


def mkdir(dir_path: Union[str, os.PathLike], exist_ok: bool = False):
    """Create a directory.

    Args:
        dir_path: The path of the directory to create.
        exist_ok: Whether to ignore the error if the directory
            already exists.
    """
    fs = get_mapper(str(dir_path)).fs
    fs.mkdirs(str(dir_path), exist_ok=exist_ok)


def md5(filepath: Union[str, os.PathLike]):
    """Return the md5 hash of a file.

    Args:
        filepath: The path to the file to compute the MD5 hash on.
    """
    with fsspec.open(filepath) as f:
        file_hash = hashlib.md5()
        file_hash.update(f.read())  # type: ignore
        file_hash = file_hash.hexdigest()
    return file_hash


def glob(path: str, detail: bool = False, **kwargs) -> List[str]:
    """Find files by glob-matching.

    Args:
        path: A glob-style path.
    """
    # Get the list of paths
    fs = get_mapper(path).fs
    paths = fs.glob(path, detail=detail, **kwargs)
    paths = [fsspec.utils._unstrip_protocol(d, fs) for d in paths]
    return paths


def copy_dir(
    source: Union[str, pathlib.Path],
    destination: Union[str, pathlib.Path],
    force: bool = False,
    progress: bool = False,
    leave_progress: bool = True,
    file_progress: bool = False,
    file_leave_progress: bool = False,
    chunk_size: Optional[int] = None,
):
    """Copy one directory to another location across different filesystem (local, S3, GCS, etc).

    Note that if both FS from source and destination are the same, progress won't be shown.

    Args:
        source: Path to the source directory.
        destination: Path to the destination directory.
        chunk_size: the chunk size to use. If progress is enabled the chunk
            size is `None`, it is set to 2048.
        force: whether to overwrite the destination directory if it exists.
        progress: Whether to display a progress bar.
        leave_progress: Whether to hide the progress bar once the copy is done.
        file_progress: Whether to display a progress bar for each file.
        file_leave_progress: Whether to hide the progress bar once a file copy is done.
        chunk_size: See `dm.utils.fs.copy_file`.
    """

    source = str(source)
    destination = str(destination)

    source_fs = get_mapper(source).fs
    destination_fs = get_mapper(destination).fs

    # Sanity check
    if not is_dir(source):
        raise ValueError(
            f"The directory being copied does not exist or is not a directory: {source}"
        )

    if not force and is_dir(destination):
        raise ValueError(f"The destination folder to copy already exists: {destination}")

    # If both fs are the same then we just rely on the internal `copy` method
    # which is much faster.
    if destination_fs.__class__ == source_fs.__class__:
        source_fs.copy(source, destination, recursive=True)
        return

    # Get all input paths with details
    # NOTE(hadim): we could have use `.glob(..., detail=True)` here but that API is inconsistent
    # between the backends resulting in different object types being returned (dict, list, etc).
    detailed_paths = source_fs.find(source, withdirs=True, detail=True)
    detailed_paths = list(detailed_paths.values())

    # Get list of input types
    input_types = [d["type"] for d in detailed_paths]

    # Get list of input path + add protocol if needed
    input_paths = [d["name"] for d in detailed_paths]
    input_paths = [fsspec.utils._unstrip_protocol(p, source_fs) for p in input_paths]

    # Build all the output paths
    output_paths: List[str] = fsspec.utils.other_paths(input_paths, destination)  # type: ignore

    def _copy_source_to_destination(input_path, input_type, output_path):
        # A directory
        if input_type == "directory":
            destination_fs.mkdir(output_path)

        # A file
        else:
            copy_file(
                input_path,
                output_path,
                force=force,
                progress=file_progress,
                leave_progress=file_leave_progress,
                chunk_size=chunk_size,
            )

    # Copy source files/directories to destination in parallel
    parallelized(
        _copy_source_to_destination,
        inputs_list=list(zip(input_paths, input_types, output_paths)),
        arg_type="args",
        progress=progress,
        tqdm_kwargs=dict(leave=leave_progress),
        scheduler="threads",
    )
