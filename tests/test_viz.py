import pytest

import numpy as np
import ipywidgets as widgets

import PIL
from PIL import Image

import datamol as dm


# NOTE(hadim): rdkit returns different image objects
# according to the Python process context (Jupyter notebook vs terminal).
# In consequence, those tests will fail if they are executed within a
# Jupyter notebook.


def test_to_image():

    # Get a list of molecules
    data = dm.data.freesolv()
    mols = dm.from_df(data)  # type: ignore
    mols = mols[:8]

    # With multiple molecules
    legends = [dm.to_smiles(mol) for mol in mols]
    image = dm.viz.to_image(mols, legends=legends, n_cols=4, mol_size=(200, 200))
    # image = _convert_ipython_to_array(image)
    image = np.array(image)

    assert image.dtype == np.uint8
    assert image.shape == (400, 800, 3)
    assert image.shape[1] == 200 * 4

    # With a single molecule
    mol = mols[0]
    legends = dm.to_smiles(mol)
    image = dm.viz.to_image(mol, legends=legends, mol_size=(200, 200))
    # image = _convert_ipython_to_array(image)
    image = np.array(image)

    assert image.dtype == np.uint8
    assert image.shape == (200, 200, 3)

    dm.viz.to_image(mol, indices=True, mol_size=400)


def test_to_image_incorrect_aromaticity():
    query = "C-c1cn(-C-2-[N,O:3]-[#6@H](-C-[#6,#8:1]-[*:2])-C(-[#8])-C-2-[#1,#8,#9:4])c2ncnc(-C)c12"
    mol = dm.from_smarts(query)
    dm.to_image(
        mol,
        mol_size=300,
        use_svg=False,
        legends="a legend",
        legend_fontsize=40,
        stereo_annotations=False,
    )


def test_to_image_save_file(tmpdir):
    smiles = "CCCOCc1cc(c2ncccc2)ccc1"
    mol = dm.to_mol(smiles)

    image_path = str(tmpdir.join("mol.png"))
    dm.viz.to_image(mol, outfile=image_path, use_svg=False)

    # check whether the png is valid
    try:
        img = Image.open(image_path)
        img.verify()
    except PIL.UnidentifiedImageError:
        pytest.fail(f"The image {image_path} is invalid.")

    image_path = str(tmpdir.join("mol.svg"))
    dm.viz.to_image(mol, outfile=image_path, use_svg=True)

    # check whether the svg looks valid
    with open(image_path) as f:
        content = f.read().strip()
    assert content.startswith("<?xml ")
    assert content.endswith("</svg>")


def test_conformers():

    import nglview as nv

    smiles = "CCCC=O"
    mol = dm.to_mol(smiles)
    mol = dm.conformers.generate(mol)

    # one conformer
    view = dm.viz.conformers(mol)
    assert type(view) == nv.widget.NGLWidget

    # multiple conformers
    view = dm.viz.conformers(mol, n_confs=12)
    assert type(view) == widgets.GridspecLayout


def test_circle_grid(tmp_path):
    mol = dm.to_mol("CC(=O)OC1=CC=CC=C1C(=O)O")

    im = dm.viz.circle_grid(
        mol,
        [
            [dm.to_mol("CCC"), dm.to_mol("CCCCCCC")],
            [dm.to_mol("CCCO"), dm.to_mol("CCCCCCCO")],
        ],
    )

    im.save(tmp_path / "image.png")


def test_match_substructure():
    mol1 = dm.to_mol("CC(=O)OC1=CC=CC=C1C(=O)O")
    mol2 = dm.to_mol("CCN(CC)CC(=O)CC(C)NC1=C2C=CC(=CC2=NC=C1)Cl")

    pattern1 = dm.from_smarts("[C;H0](=O)")
    pattern2 = dm.to_mol("CN(C)")

    # Test multiple scenarios

    dm.viz.match_substructure(
        mols=[mol1, mol2],
        patterns=[pattern1, pattern2],
        highlight_bonds=True,
        use_svg=True,
    )
    dm.viz.match_substructure(
        mols=mol1,
        patterns=[pattern1, pattern2],
        highlight_bonds=True,
        use_svg=True,
    )
    dm.viz.match_substructure(
        mols=[mol1, mol2],
        patterns=pattern1,
        highlight_bonds=False,
        use_svg=False,
    )
    dm.viz.match_substructure(
        mols=mol1,
        patterns=pattern2,
        highlight_bonds=True,
        use_svg=False,
    )
