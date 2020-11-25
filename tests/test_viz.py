import unittest

import numpy as np
import ipywidgets as widgets

import datamol as dm


class TestViz(unittest.TestCase):
    def test_to_image(self):

        # Get a list of molecules
        data = dm.data.freesolv()
        mols = dm.from_df(data)
        mols = mols[:8]

        # With multiple molecules
        legends = [dm.to_smiles(mol) for mol in mols]
        image = dm.viz.to_image(mols, legends=legends, n_cols=4, mol_size=(200, 200))

        image = np.array(image)
        assert image.dtype == np.uint8
        assert image.shape == (400, 800, 3)
        assert image.shape[1] == 200 * 4

        # With a single molecule
        mol = mols[0]
        legends = dm.to_smiles(mol)
        image = dm.viz.to_image(mol, legends=legends, mol_size=(200, 200))

        image = np.array(image)
        print(image.shape)
        assert image.dtype == np.uint8
        assert image.shape == (200, 200, 3)

    def test_conformers(self):

        import nglview as nv

        smiles = "CCCC=O"
        mol = dm.to_mol(smiles)
        mol = dm.conformers.generate(mol, n_confs=2)

        # one conformer
        view = dm.viz.conformers(mol)
        assert type(view) == nv.widget.NGLWidget

        # multiple conformers
        view = dm.viz.conformers(mol, n_confs=12)
        assert type(view) == widgets.GridspecLayout
