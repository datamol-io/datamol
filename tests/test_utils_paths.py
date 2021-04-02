import pathlib
import tempfile

import datamol as dm


def test_copy_files():
    source_path = pathlib.Path(tempfile.NamedTemporaryFile(delete=False).name)
    destination_path = pathlib.Path(tempfile.NamedTemporaryFile(delete=False).name)

    content = "hello this is a content"
    with open(source_path, "w") as f:
        f.write(content)

    dm.utils.fs.copy_file(source_path, destination_path)

    with open(destination_path) as f:
        f.read() == content

    source_path.unlink()
    destination_path.unlink()
