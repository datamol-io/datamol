import pytest
import os
import pandas as pd
from pandas.testing import assert_frame_equal
import datamol as dm

@pytest.fixture
def tmp_path(tmpdir):
    return str(tmpdir.mkdir("temp").join("test_file"))
    
def test_dataframe_csv(tmp_path):
    data = pd.DataFrame({"col1": [1, 2, 3], "col2": ["a", "b", "c"]})
    dm.save_dataframe(data, tmp_path + ".csv")
    loaded_data = dm.open_dataframe(tmp_path + ".csv")
    assert_frame_equal(data, loaded_data)

def test_dataframe_excel(tmp_path):
    data = pd.DataFrame({"col1": [1, 2, 3], "col2": ["a", "b", "c"]})
    dm.save_dataframe(data, tmp_path + ".xlsx")
    loaded_data = dm.open_dataframe(tmp_path + ".xlsx")
    assert_frame_equal(data, loaded_data)

def test_dataframe_parquet(tmp_path):
    data = pd.DataFrame({"col1": [1, 2, 3], "col2": ["a", "b", "c"]})
    dm.save_dataframe(data, tmp_path + ".parquet")
    loaded_data = dm.open_dataframe(tmp_path + ".parquet")
    assert_frame_equal(data, loaded_data)

def test_dataframe_json(tmp_path):
    data = pd.DataFrame({"col1": [1, 2, 3], "col2": ["a", "b", "c"]})
    dm.save_dataframe(data, tmp_path + ".json")
    loaded_data = dm.open_dataframe(tmp_path + ".json")
    assert_frame_equal(data, loaded_data)

def test_dataframe_sdf(tmp_path):
    data = pd.DataFrame({"smiles": ["CC", "CCCC", "CCC"]})
    dm.save_dataframe(data, tmp_path + ".sdf")
    loaded_data = dm.open_dataframe(tmp_path + ".sdf")
    assert_frame_equal(data, loaded_data)

def test_save_dataframe_invalid(tmp_path):
    data = pd.DataFrame({"col1": [1, 2, 3], "col2": ["a", "b", "c"]})
    with pytest.raises(ValueError):
        dm.save_dataframe(data, tmp_path + ".invalid")
        
def test_load_dataframe_invalid(tmp_path):
    with pytest.raises(ValueError):
        dm.open_dataframe(tmp_path + ".invalid")
