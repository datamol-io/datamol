from typing import cast

import pandas as pd
import datamol as dm

EXTENSIONS_DICT = {
    "csv": [
        ".csv",
        ".csv.gz",
        ".csv.bz2",
        ".csv.zip",
        ".csv.xz",
        ".csv.zst",
        ".csv.tar",
        ".csv.tar.gz",
        ".csv.tar.xz",
        ".csv.tar.bz2",
    ],
    "excel": [".xlsx"],
    "parquet": [".parquet"],
    "json": [
        ".json",
        ".json.gz",
        ".json.bz2",
        ".json.zip",
        ".json.xz",
        ".json.zst",
        ".json.tar",
        ".json.tar.gz",
        ".json.tar.xz",
        ".json.tar.bz2",
    ],
    "sdf": [".sdf"],
}


def _guess_filetype(path: str):
    """Return a filetype given an input path. Filetypes returned can be from
    `csv, excel, parquet, json, sdf`.
    """
    for name, extensions in EXTENSIONS_DICT.items():
        for ext in extensions:
            if path.endswith(ext):
                return name


def open_dataframe(path: str) -> pd.DataFrame:
    """Open a dataframe file whatever its filetype from
    `csv, excel, parquet, json, sdf`.
    
    args:
        path: path to the file.
    """

    filetype = _guess_filetype(path)

    data = None
    if filetype == "csv":
        data = pd.read_csv(path)
    elif filetype == "excel":
        data = pd.read_excel(path)
    elif filetype == "parquet":
        data = pd.read_parquet(path)
    elif filetype == "json":
        data = pd.read_json(path)
    elif filetype == "sdf":
        data = dm.read_sdf(path, as_df=True)

    if data is None:
        raise ValueError(f"The file type of {path} is not supported.")

    data = cast(pd.DataFrame, data)

    return data


def save_dataframe(data: pd.DataFrame, path: str):
    """Save a dataframe file whatever its filetype from
    `csv, excel, parquet, json, sdf`.
    
    args:
        data: dataframe to save.
        path: path to save the file.
    """

    filetype = _guess_filetype(path)

    if filetype == "csv":
        data.to_csv(path, index=False)
    elif filetype == "excel":
        data.to_excel(path, index=False)
    elif filetype == "parquet":
        data.to_parquet(path)
    elif filetype == "json":
        data.to_json(path)
    elif filetype == "sdf":
        dm.to_sdf(data, path)
    else:
        raise ValueError(f"The file type of {path} is not supported.")
