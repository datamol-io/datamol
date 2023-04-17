import pytest

import datamol as dm


def test_datamol_import_fails():
    with pytest.raises(AttributeError):
        dm.that_import_does_not_exist
