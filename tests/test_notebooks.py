import pytest
import pathlib

import nbformat
from nbconvert.preprocessors import ExecutePreprocessor

ROOT_DIR = pathlib.Path(__file__).parent.resolve()

NOTEBOOK_DIR = ROOT_DIR.parent / "docs" / "tutorials"

NOTEBOOK_PATHS = NOTEBOOK_DIR.glob("*.ipynb")
NOTEBOOK_PATHS = sorted(list(NOTEBOOK_DIR.glob("*.ipynb")))


@pytest.mark.skip_platform("win")
@pytest.mark.parametrize("nb_path", NOTEBOOK_PATHS, ids=[str(n.name) for n in NOTEBOOK_PATHS])
def test_notebook(nb_path):

    # Setup and configure the processor to execute the notebook
    ep = ExecutePreprocessor(timeout=600, kernel_name="python")

    # Open the notebook
    with open(nb_path) as f:
        nb = nbformat.read(f, as_version=nbformat.NO_CONVERT)

    # Execute the notebook
    ep.preprocess(nb, {"metadata": {"path": NOTEBOOK_DIR}})
