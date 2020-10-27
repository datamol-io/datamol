import pathlib
import shutil

from loguru import logger

import types
import pathlib
import importlib

import pandas as pd
import datamol


def _introspect(module, parent_lib):

    data = []
    for attr_str in dir(module):

        if attr_str.startswith("_"):
            continue

        # Get the attribute
        attr = getattr(module, attr_str)

        # Attribute type
        attr_type = type(attr)
        attr_type_str = attr_type.__name__

        # Get the attribute's parent module (the first level module after `datamol`)
        if isinstance(attr, types.ModuleType):
            parent_module_str = [module.__name__]
            is_module = True
        else:
            parent_module_str = attr.__module__.split(".")
            is_module = False

        datum = {}
        datum["attr"] = attr_str
        datum["attr_type"] = attr_type_str
        datum["parent_list"] = parent_module_str
        datum["parent"] = parent_module_str[-1]
        datum["is_module"] = is_module
        data.append(datum)

    # Build dataframe
    data = pd.DataFrame(data)
    return data


def datamol_introspector(module, parent_lib):
    """This function is a complex and VERY opinated way to
    inspect a given Python module in order to generate
    an API documentation from it.
    """

    data = _introspect(module, parent_lib=parent_lib)

    # Filter out attribute outside of the parent library
    mask = data["parent_list"].apply(lambda x: x[0] == parent_lib)
    data = data[mask]

    # Sanity check
    assert data["attr_type"].unique().tolist() == ["type", "module", "function"]

    # Get attributes by their types
    mask = data["attr_type"] == "module"
    all_non_modules = data[~mask]
    all_modules = data[mask]

    # Now we remove from the modules the one already in `all_non_modules`
    to_remove = all_non_modules["parent"].unique().tolist()
    mask = all_modules["attr"].isin(to_remove)
    all_modules = all_modules[~mask]

    return all_modules, all_non_modules


def generate_mkdocs_api(docs_dir):
    """Create an `api` folder in `docs` (the MkDocs `docs_dir` config).
    In this folder, generate the API markdown files of datamol.
    """

    # Prepare the API doc directory
    docs_dir = pathlib.Path(docs_dir)
    api_dir = docs_dir / "api"
    logger.info(f"api_dir: {api_dir}")

    if api_dir.is_dir():
        shutil.rmtree(api_dir)
    api_dir.mkdir()

    # Get attributes data
    parent_lib = "datamol"
    all_modules, all_non_modules = datamol_introspector(datamol, parent_lib="datamol")

    # `datamol.md`
    with open(api_dir / f"{parent_lib}.md", "w") as f:
        f.write(f"# `{parent_lib}`\n\n")

        f.write(f"All the below functions are accessible under `datamol.FUNCTION_NAME`.\n\n")

        for _, rows in all_non_modules.groupby("parent"):
            parent = rows["parent"].iloc[0]
            f.write(f"## `{parent_lib}.{parent}`\n\n")

            for _, attr in rows.iterrows():
                f.write(f"::: {parent_lib}:{attr['attr']}\n")
            f.write("\n")

    # `datamol.SUBMODULE.md`
    for _, row in all_modules.iterrows():
        with open(api_dir / f"{row['parent']}.{row['attr']}.md", "w") as f:
            f.write(f"# `{row['parent']}.{row['attr']}`\n\n")

            f.write(
                f"All the below functions are accessible under `{row['parent']}.{row['attr']}`.\n\n"
            )

            # Get data from this module
            module = importlib.import_module(f"{row['parent']}.{row['attr']}")
            data = _introspect(module, parent_lib)

            # Sanity check
            #         assert data["attr_type"].unique().tolist() == ['type', 'function']

            for _, attr_child in data.iterrows():
                f.write(f"::: {row['parent']}.{row['attr']}:{attr_child['attr']}\n")
            f.write("\n")

    logger.info(f"API doc generation is done.")


if __name__ == "__main__":
    generate_mkdocs_api("./docs/")
