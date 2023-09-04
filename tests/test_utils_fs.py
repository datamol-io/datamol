import pytest

import pathlib

import fsspec
import datamol as dm


def test_copy_files(tmp_path):
    source_path = tmp_path / "source.txt"
    destination_path = tmp_path / "destination.txt"

    content = "hello this is a content"
    with open(source_path, "w") as f:
        f.write(content)

    dm.utils.fs.copy_file(source_path, destination_path)

    with open(destination_path) as f:
        f.read() == content


def test_copy_dir(tmp_path):
    source_path = tmp_path / "source_dir"
    source_path_subdir = source_path / "a_subdir"
    destination_path = tmp_path / "destination_dir"
    destination_path_subdir = destination_path / "a_subdir"

    dm.utils.fs.mkdir(source_path)
    dm.utils.fs.mkdir(source_path_subdir)

    content = "hello this is a content"
    file1_path = source_path / "hello.txt"
    with open(file1_path, "w") as f:
        f.write(content)

    file2_path = source_path_subdir / "hello.txt"
    with open(file2_path, "w") as f:
        f.write(content)

    assert not dm.utils.fs.is_dir(destination_path_subdir)
    assert not dm.utils.fs.is_dir(destination_path)

    dm.utils.fs.copy_dir(source_path, destination_path)

    assert dm.utils.fs.is_dir(destination_path_subdir)
    assert dm.utils.fs.is_dir(destination_path)
    assert dm.utils.fs.is_file(file1_path)
    assert dm.utils.fs.is_file(file2_path)

    with open(file1_path) as f:
        f.read() == content

    with open(file2_path) as f:
        f.read() == content


def test_mkdir(tmp_path):
    source_path = tmp_path / "source_dir"
    source_path_subdir = source_path / "a_subdir"

    dm.utils.fs.mkdir(source_path)

    assert dm.utils.fs.is_dir(source_path)
    assert not dm.utils.fs.is_dir(source_path_subdir)

    dm.utils.fs.mkdir(source_path_subdir)

    assert dm.utils.fs.is_dir(source_path)
    assert dm.utils.fs.is_dir(source_path_subdir)


@pytest.mark.skip_platform("win")
def test_cache_dir():
    cache_dir = dm.utils.fs.get_cache_dir("my_app")
    assert str(cache_dir).endswith("my_app")
    assert cache_dir.exists()
    assert cache_dir.is_dir()

    cache_dir = dm.utils.fs.get_cache_dir("my_app", suffix="likelydonotalreadyexist", create=False)
    assert str(cache_dir).endswith("likelydonotalreadyexist")
    assert not cache_dir.exists()
    assert not cache_dir.is_dir()

    cache_dir = dm.utils.fs.get_cache_dir("my_app", suffix="iamasuffix")
    assert str(cache_dir).endswith("iamasuffix")
    assert "my_app" in str(cache_dir)
    assert cache_dir.exists()
    assert cache_dir.is_dir()


def test_get_mapper(tmp_path):
    fsmapper = dm.utils.fs.get_mapper(str(tmp_path / "test.txt"))
    assert fsmapper.fs.protocol == "file"


@pytest.mark.skip_platform("win")
def test_get_basename(tmp_path):
    assert dm.utils.fs.get_basename(str(tmp_path / "test.txt")) == "test.txt"
    assert dm.utils.fs.get_basename("s3://a-bucket-that-likely-do-not-exist/test.txt") == "test.txt"


def test_get_extension(tmp_path):
    assert dm.utils.fs.get_extension(str(tmp_path / "test.txt")) == "txt"
    assert dm.utils.fs.get_extension("s3://a-bucket-that-likely-do-not-exist/test.txt") == "txt"


def test_exists(tmp_path):
    tmp_file = tmp_path / "test.txt"

    assert not dm.utils.fs.exists(tmp_file)
    assert not dm.utils.fs.is_file(tmp_file)

    assert dm.utils.fs.is_dir(tmp_path)
    assert not dm.utils.fs.is_dir(tmp_path / "likely-does-not-exist")

    with open(tmp_file, "w") as f:
        f.write("hello")

    assert dm.utils.fs.exists(tmp_file)
    assert dm.utils.fs.is_file(tmp_file)

    assert not dm.utils.fs.is_file(open(tmp_file))
    assert not dm.utils.fs.is_dir(open(tmp_file))


def test_get_protocol(tmp_path):
    assert dm.utils.fs.get_protocol(tmp_path / "ahahah.txt") == "file"
    assert dm.utils.fs.get_protocol("s3://a-bucket-that-likely-do-not-exist/test.txt") == "s3"


def test_is_local_path(tmp_path):
    assert dm.utils.fs.is_local_path(tmp_path / "ahahah.txt")
    assert not dm.utils.fs.is_local_path("s3://a-bucket-that-likely-do-not-exist/test.txt")


@pytest.mark.skip_platform("win")
def test_join(tmp_path):
    assert (
        dm.utils.fs.join("s3://a-bucket-that-likely-do-not-exist", "test.txt")
        == "s3://a-bucket-that-likely-do-not-exist/test.txt"
    )
    assert dm.utils.fs.join(tmp_path, "test.txt") == str(tmp_path / "test.txt")


def test_get_size(tmp_path):
    tmp_file = tmp_path / "test.txt"

    with open(tmp_file, "w") as f:
        f.write("hello")

    assert dm.utils.fs.get_size(tmp_file) > 0
    assert dm.utils.fs.get_size(open(tmp_file)) > 0
    assert dm.utils.fs.get_size(fsspec.open(tmp_file)) > 0


def test_md5(tmp_path):
    tmp_file = tmp_path / "test.txt"

    with open(tmp_file, "w") as f:
        f.write("hello")

    assert dm.utils.fs.md5(tmp_file) == "5d41402abc4b2a76b9719d911017c592"


@pytest.mark.skip_platform("win")
def test_glob(tmp_path):
    for i in range(5):
        tmp_file = tmp_path / f"test_{i}.txt"

        with open(tmp_file, "w") as f:
            f.write("hello")

    tmp_path_regex = tmp_path / "*.txt"
    assert len(dm.utils.fs.glob(tmp_path_regex)) == 5


def test_copy_file(tmp_path):
    tmp_file = tmp_path / "test.txt"

    assert dm.utils.fs.is_dir(tmp_path)
    assert dm.utils.fs.is_dir(str(tmp_path))
    assert dm.utils.fs.is_dir(pathlib.Path(str(tmp_path)))

    assert not dm.utils.fs.is_dir(tmp_path / "not_exist_dir")
    assert not dm.utils.fs.is_dir(str(tmp_path / "not_exist_dir"))
    assert not dm.utils.fs.is_dir(pathlib.Path(str(tmp_path / "not_exist_dir")))

    with open(tmp_file, "w") as f:
        f.write("hello")

    tmp_file2 = tmp_path / "test2.txt"
    assert not dm.utils.fs.is_file(tmp_file2)
    assert not dm.utils.fs.is_file(str(tmp_file2))
    assert not dm.utils.fs.is_file(pathlib.Path(str(tmp_file2)))

    dm.utils.fs.copy_file(tmp_file, tmp_file2)

    assert dm.utils.fs.is_file(tmp_file2)
    assert dm.utils.fs.is_file(str(tmp_file2))
    assert dm.utils.fs.is_file(pathlib.Path(str(tmp_file2)))
    assert open(tmp_file2).read() == "hello"

    with pytest.raises(ValueError):
        dm.utils.fs.copy_file(tmp_file, tmp_file2)

    tmp_file3 = tmp_path / "test3.txt"
    dm.utils.fs.copy_file(tmp_file, tmp_file3, progress=True)
    assert dm.utils.fs.is_file(tmp_file3)
    assert dm.utils.fs.is_file(str(tmp_file3))
    assert dm.utils.fs.is_file(pathlib.Path(str(tmp_file3)))
    assert open(tmp_file3).read() == "hello"
