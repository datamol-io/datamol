import pytest

import pandas as pd
import datamol as dm


def test_template_align():
    data: pd.DataFrame = dm.cdk2(as_df=True)  # type: ignore
    data = data.iloc[:6].copy()  # type: ignore

    template = data.iloc[0]["mol"]
    data["aligned_mol"] = data["mol"].apply(lambda x: dm.align.template_align(x, template=template))
    assert bool(data["aligned_mol"].apply(lambda x: isinstance(x, dm.Mol)).all()) is True

    template = data.iloc[0]["smiles"]
    data["aligned_mol"] = data["smiles"].apply(
        lambda x: dm.align.template_align(x, template=template)
    )
    assert bool(data["aligned_mol"].apply(lambda x: isinstance(x, dm.Mol)).all()) is True

    template = data.iloc[0]["mol"]
    data["aligned_mol"] = data["mol"].apply(
        lambda x: dm.align.template_align(x, template=template, auto_select_coord_gen=True)
    )
    assert bool(data["aligned_mol"].apply(lambda x: isinstance(x, dm.Mol)).all()) is True

    template = data.iloc[0]["mol"]
    data["aligned_mol"] = data["mol"].apply(
        lambda x: dm.align.template_align(x, template=template, use_depiction=False)
    )
    assert bool(data["aligned_mol"].apply(lambda x: isinstance(x, dm.Mol)).all()) is True

    template = None
    data["aligned_mol"] = data["mol"].apply(lambda x: dm.align.template_align(x, template=template))
    assert bool(data["aligned_mol"].apply(lambda x: isinstance(x, dm.Mol)).all()) is True

    template = None
    data["aligned_mol"] = data["mol"].apply(
        lambda x: dm.align.template_align(x, template=template, copy=False)
    )
    assert bool(data["aligned_mol"].apply(lambda x: isinstance(x, dm.Mol)).all()) is True

    assert dm.align.template_align(None) is None


def test_auto_align_many():
    data: pd.DataFrame = dm.solubility(as_df=True)  # type: ignore
    data = data.iloc[:16].copy()  # type: ignore

    excepted_cluster_size = [8, 6, 5, 6, 6]

    for i, partition_method in enumerate(
        [
            "cluster",
            "scaffold",
            "anongraph-scaffold",
            "anon-scaffold",
            "strip-scaffold",
        ]
    ):
        print(partition_method)

        data["aligned_mol"] = dm.align.auto_align_many(
            data["mol"],
            partition_method=partition_method,
        )

        props = data["aligned_mol"].apply(lambda x: pd.Series(x.GetPropsAsDict()))

        assert "dm.auto_align_many.cluster_id" in props.columns
        assert "dm.auto_align_many.core" in props.columns
        assert props["dm.auto_align_many.cluster_id"].dtype.name == "int64"
        assert props["dm.auto_align_many.core"].dtype.name == "object"

        assert props["dm.auto_align_many.cluster_id"].unique().shape[0] == excepted_cluster_size[i]

    with pytest.raises(ValueError):
        dm.align.auto_align_many(data["mol"], partition_method="invalid")
