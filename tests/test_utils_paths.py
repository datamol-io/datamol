import pathlib
import tempfile

import datamol as dm


def test_copy_files():
    source_path = pathlib.Path(tempfile.mkstemp()[1])
    destination_path = pathlib.Path(tempfile.mkstemp()[1])

    content = "hello this is a content"
    with open(source_path, "w") as f:
        f.write(content)

    dm.utils.copy_files(source_path, destination_path)

    with open(destination_path) as f:
        f.read() == content

    source_path.unlink()
    destination_path.unlink()
