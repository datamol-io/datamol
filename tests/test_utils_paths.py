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
