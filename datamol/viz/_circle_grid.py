from typing import List
from typing import Union
from typing import Optional
from typing import Any

import math
import numpy as np
import random
from rdkit.Chem import Draw
from rdkit.Geometry.rdGeometry import Point2D

from .utils import drawer_to_image
from .utils import prepare_mol_for_drawing
from .utils import image_to_file
from .utils import to_rdkit_color

from datamol.types import DatamolColor
from datamol.types import Mol
from datamol.types import RDKitColor

import datamol as dm


def circle_grid(
    center_mol: Union[Mol, str],
    ring_mols: List[List[Union[Mol, str]]],
    act_mapper: Optional[dict] = None,
    margin: int = 50,
    legend: Optional[str] = None,
    ring_scaler: float = 1.0,
    align: Optional[Union[Mol, str, bool]] = None,
    use_svg: bool = True,
    ring_color: Optional[DatamolColor] = None,
    ring_mol_start_angles_degrees: Optional[List[float]] = None,
    center_mol_highlight_atoms: Optional[List[int]] = None,
    center_mol_highlight_bonds: Optional[List[int]] = None,
    ring_mol_highlight_atoms: Optional[List[List[int]]] = None,
    ring_mol_highlight_bonds: Optional[List[List[int]]] = None,
    outfile: Optional[str] = None,
    kekulize: bool = True,
    layout_random_seed: Optional[int] = 19,
    **kwargs: Any,
):
    """Show molecules in concentric rings, with one molecule at the center
    Args:
        center_mol: Molecule at center of the rings
        ring_mols: List of molecule for each level of concentric rings around the center mol
        legend: optional global legend for the figure
        margin: Margin between the circle layers
        ring_scaler: Scale the size of the molecules in each circle by this factor compared to the center molecule
        act_mapper: dictionary of activity for each molecule
        align: Whether to align the 2D coordinates of the molecules.
            - If set to True, align all molecules with `dm.align.auto_align_many()`.
            - If set to a molecule, it is used as a template for alignment with `dm.align.template_align()`.
            - If set to False, no alignment is performed.
            - If set to None (default), the ring (peripheral) molecules are aligned to the center molecules.
        use_svg: Whether to use SVG or use PNG
        center_mol_highlight_atoms: List of atom indices to highlight for the center molecule
        center_mol_highlight_bonds: List of bond indices to highlight for the center molecule
        ring_mol_highlight_atoms: List of list of atom indices to highlight for molecules at each level of the concentric rings
        ring_mol_highlight_bonds: List of list of bond indices to highlight for molecules at each level of the concentric rings
        ring_color: Color of the concentric rings. Set to None to not draw any ring.
        ring_mol_start_angles_degrees: List of angles in degrees to start drawing the molecules at each level of the concentric
            rings. If None then a random position will be used.
        kekulize: Whether to kekulize the molecules before drawing.
        layout_random_seed: Random seed for the layout of the molecules. Set to None for no seed.
        outfile: Optional path to the save the output file.
        **kwargs: Additional arguments to pass to the drawing function. See RDKit
            documentation related to `MolDrawOptions` for more details at
            https://www.rdkit.org/docs/source/rdkit.Chem.Draw.rdMolDraw2D.html.

    """

    grid = MolsCircleGrid(
        center_mol,
        ring_mols,
        legend=legend,
        margin=margin,
        ring_scaler=ring_scaler,
        act_mapper=act_mapper,
        align=align,
        use_svg=use_svg,
        ring_color=ring_color,
        ring_mol_start_angles_degrees=ring_mol_start_angles_degrees,
        center_mol_highlight_atoms=center_mol_highlight_atoms,
        center_mol_highlight_bonds=center_mol_highlight_bonds,
        ring_mol_highlight_atoms=ring_mol_highlight_atoms,
        ring_mol_highlight_bonds=ring_mol_highlight_bonds,
        kekulize=kekulize,
        layout_random_seed=layout_random_seed,
        **kwargs,
    )
    return grid(outfile=outfile)


class MolsCircleGrid:
    def __init__(
        self,
        center_mol: Union[Mol, str],
        ring_mols: List[List[Union[Mol, str]]],
        legend: Optional[str] = None,
        margin: int = 50,
        ring_scaler: float = 1.0,
        act_mapper: Optional[dict] = None,
        align: Optional[Union[Mol, str, bool]] = None,
        use_svg: bool = True,
        line_width: Optional[float] = None,
        ring_color: Optional[DatamolColor] = None,
        ring_mol_start_angles_degrees: Optional[List[float]] = None,
        center_mol_highlight_atoms: Optional[List[int]] = None,
        center_mol_highlight_bonds: Optional[List[int]] = None,
        ring_mol_highlight_atoms: Optional[List[List[int]]] = None,
        ring_mol_highlight_bonds: Optional[List[List[int]]] = None,
        kekulize: bool = True,
        layout_random_seed: Optional[int] = 19,
        **kwargs: Any,
    ):
        """Show molecules in concentric rings, with one molecule at the center

        Args:
            center_mol: Molecule at center of the rings
            ring_mols: List of molecule for each level of concentric rings around the center mol
            legend: optional global legend for the figure
            margin: Margin between the circle layers
            ring_scaler: Scale the size of the molecules in each circle by this factor compared to the center molecule
            act_mapper: dictionary of activity for each molecule
            align: Whether to align the 2D coordinates of the molecules.
                - If set to True, align all molecules with `dm.align.auto_align_many()`.
                - If set to a molecule, it is used as a template for alignment with `dm.align.template_align()`.
                - If set to False, no alignment is performed.
                - If set to None (default), the ring (peripheral) molecules are aligned to the center molecules. This is the default behaviour.
            use_svg: Whether to use SVG or use PNG
            line_width: Width of the lines to draw
            center_mol_highlight_atoms: List of atom indices to highlight for the center molecule
            center_mol_highlight_bonds: List of bond indices to highlight for the center molecule
            ring_mol_highlight_atoms: List of list of atom indices to highlight for molecules at each level of the concentric rings
            ring_mol_highlight_bonds: List of list of bond indices to highlight for molecules at each level of the concentric rings
            ring_color: Color of the concentric rings. Set to None to not draw any ring.
            ring_mol_start_angles_degrees: List of angles in degrees to start drawing the molecules at each level of the concentric
                rings. If None then a random position will be used.
            kekulize: Whether to kekulize the molecules before drawing.
            layout_random_seed: Random seed for the layout of the molecules. Set to None for no seed.
            **kwargs: Additional arguments to pass to the drawing function. See RDKit
                documentation related to `MolDrawOptions` for more details at
                https://www.rdkit.org/docs/source/rdkit.Chem.Draw.rdMolDraw2D.html.
        """
        assert dm.is_greater_than_current_rdkit_version(
            "2023.03"
        ), "MolsCircleGrid requires RDKit version >= 2022.09"
        self.center_mol = center_mol
        self.ring_mols = ring_mols
        if not isinstance(self.ring_mols[0], (tuple, list, np.ndarray)):
            self.ring_mols = [self.ring_mols]
        self.ring_count = len(self.ring_mols)
        self.legend = legend or None
        self.margin = margin
        self.ring_scaler = ring_scaler
        self.align = align
        self.use_svg = use_svg
        self.line_width = line_width
        self.ring_color = ring_color
        self.ring_mol_start_angles_degrees = ring_mol_start_angles_degrees
        self.ring_color_rdkit: Optional[RDKitColor] = to_rdkit_color(ring_color)
        self.ring_mol_highlight_atoms = ring_mol_highlight_atoms
        self.ring_mol_highlight_bonds = ring_mol_highlight_bonds
        self.center_mol_highlight_atoms = center_mol_highlight_atoms
        self.center_mol_highlight_bonds = center_mol_highlight_bonds
        self.kekulize = kekulize
        self.layout_random_seed = layout_random_seed
        self._global_legend_size = 0
        if self.legend is not None:
            self._global_legend_size = max(25, self.margin)
        self._format_activity_mapper(act_mapper)
        self._initialize_drawing(**kwargs)

    def _format_activity_mapper(self, act_mapper: Optional[dict] = None) -> None:
        """Format the activity mapper to a dictionary of legends

        Args:
            act_mapper: dictionary of activity for each molecule
        """
        self._legends = {}
        self.act_mapper = {}
        if act_mapper is None:
            return
        m_0 = list(act_mapper.keys())[0]
        if isinstance(m_0, Mol) or dm.to_mol(m_0) is not None:
            act_mapper = {dm.unique_id(k): v for k, v in act_mapper.items()}
        self.act_mapper = act_mapper
        for k, v in self.act_mapper.items():
            txt = []
            for prop, propval in v.items():
                if not isinstance(propval, str):
                    propval = "{:.2f}".format(propval)
                txt.append(f"{prop}: {propval}")
            txt = "\n".join(txt)
            self._legends[k] = txt

    def _initialize_drawing(self, **drawing_options) -> None:
        """Initialize the canvas for drawing and perform all preparation for input molecules

        Args:
            drawing_options: drawing options for the canvas
        """

        # prepare and align all the molecules first
        all_mols = [self.center_mol] + [mol for mols in self.ring_mols for mol in mols]
        all_mols = [dm.to_mol(mol) if isinstance(mol, str) else mol for mol in all_mols]
        if self.align is None:
            # align to center mol
            all_mols = [
                dm.align.template_align(mol, template=self.center_mol, use_depiction=False)
                for mol in all_mols
            ]
        elif isinstance(self.align, (dm.Mol, str)):
            # align to template
            all_mols = [
                dm.align.template_align(mol, template=self.align, use_depiction=False)
                for mol in all_mols
            ]
        elif self.align is True:
            all_mols = dm.align.auto_align_many(all_mols)

        # prepare all the molecules for drawing
        all_mols = [prepare_mol_for_drawing(mol, kekulize=self.kekulize) for mol in all_mols]

        # put back center and peripheral mols after preparing them for drawing
        self.center_mol = all_mols[0]
        cur_circle_start = 1
        for i, mols in enumerate(self.ring_mols):
            self.ring_mols[i] = all_mols[cur_circle_start : cur_circle_start + len(mols)]
            cur_circle_start += len(mols)

        # get the size of each molecule, including the scaling of the ring molecules
        self._size_info = {}
        sizer = Draw.rdMolDraw2D.MolDraw2DSVG(-1, -1)
        sizer_options = sizer.drawOptions()
        current_scaling_factor = sizer_options.scalingFactor
        for ind, mol in enumerate(all_mols):
            if ind > 0:
                sizer_options.scalingFactor = current_scaling_factor * self.ring_scaler
            mol_id = dm.unique_id(mol)
            self._size_info[mol_id] = sizer.GetMolSize(mol, legend=self._legends.get(mol_id, ""))

        # compute the minimum size of the canvas based on the size of the molecules
        max_width = max([x[0] for x in self._size_info.values()])
        max_height = max([x[1] for x in self._size_info.values()])
        self.mol_size = (max_width, max_height)
        self.size = (
            max_width + (self.margin + max_width) * 2 * self.ring_count,
            max_height + (self.margin + max_height) * 2 * self.ring_count,
        )
        self.size = (
            self.size[0],
            self.size[1] + self._global_legend_size * 2,
        )
        self.midpoint = Point2D(self.size[0] // 2, self.size[1] // 2)

        if self.use_svg:
            self.canvas = Draw.rdMolDraw2D.MolDraw2DSVG(*self.size)
        else:
            self.canvas = Draw.rdMolDraw2D.MolDraw2DCairo(*self.size)
        # use fixed parameters for bond length and font size.
        self.canvas.SetFlexiMode(True)
        self.canvas.ClearDrawing()

        # Setting the drawing options
        self.draw_options = self.canvas.drawOptions()
        for k, v in drawing_options.items():
            if hasattr(self.draw_options, k):
                setattr(self.draw_options, k, v)
            else:
                raise ValueError(f"Unknown drawing option: {k}")
        if self.line_width is not None:
            self.canvas.SetLineWidth(self.line_width)
        self.canvas.SetFillPolys(False)

    def __call__(self, outfile: Optional[str] = None):
        """Draw the circular molecule and save to file if outfile is provided

        Args:
            outfile: output file name
        """
        image = self.draw()
        if outfile is not None:
            image_to_file(image, outfile, as_svg=self.use_svg)
        return image

    def draw(self):
        """Create and draw the circular molecule image"""
        radius_list = self._draw_circles()

        # draw the center mol
        self._draw_mol_at(
            self.center_mol,
            self.midpoint,
            highlight_atom=self.center_mol_highlight_atoms,
            highlight_bond=self.center_mol_highlight_bonds,
        )

        rng = random.Random(self.layout_random_seed)

        # draw the ring mols
        self.draw_options.scalingFactor *= self.ring_scaler
        for i, mols in enumerate(self.ring_mols):
            radius = radius_list[i]
            ni = len(mols)

            if self.ring_mol_start_angles_degrees is not None:
                rand_unit = np.deg2rad(self.ring_mol_start_angles_degrees[i])
            else:
                rand_unit = rng.random() * 2 * math.pi

            for k, mol in enumerate(mols):
                center_x = radius * math.cos(2 * k * math.pi / ni + rand_unit) + self.midpoint.x
                center_y = radius * math.sin(2 * k * math.pi / ni + rand_unit) + self.midpoint.y
                center = Point2D(center_x, center_y)
                ring_atom_highlight = None
                if self.ring_mol_highlight_atoms is not None:
                    ring_atom_highlight = self.ring_mol_highlight_atoms[i][k]
                ring_bond_highlight = None
                if self.ring_mol_highlight_bonds is not None:
                    ring_bond_highlight = self.ring_mol_highlight_bonds[i][k]
                self._draw_mol_at(
                    mol,
                    center,
                    highlight_atom=ring_atom_highlight,
                    highlight_bond=ring_bond_highlight,
                )

        self.draw_options.scalingFactor /= self.ring_scaler
        # draw global legend if there is any
        if self.legend is not None:
            self.canvas.SetFontSize(self.canvas.FontSize() * 3)
            text_position = Point2D(
                self.size[0] // 2, self.size[-1] - self._global_legend_size // 2
            )
            self.canvas.DrawString(self.legend, text_position, align=0, rawCoords=True)

        self.canvas.FinishDrawing()
        return drawer_to_image(self.canvas)

    def _draw_circles(self):
        """Drawing the circular rings around the center molecule"""
        if self.ring_count <= 0:
            return []
        radius_step = (min(self.size) - min(self.mol_size) - 2 * self._global_legend_size) // (
            self.ring_count * 2
        )
        radius_list = []
        full_range = range(0, min(self.size) // 2, radius_step)
        for _, radius in enumerate(full_range):
            radius += self.margin // 2
            if radius > self.margin:
                if self.ring_color_rdkit is not None:
                    self.canvas.SetColour(self.ring_color_rdkit)
                    self.canvas.DrawArc(self.midpoint, radius, 0, 360, rawCoords=True)
            radius_list.append(radius + radius_step)
        return radius_list

    def _draw_mol_at(
        self,
        mol: dm.Mol,
        mol_center: Point2D,
        highlight_atom: Optional[List[int]] = None,
        highlight_bond: Optional[List[int]] = None,
    ):
        """Draw molecule at a given position

        Args:
            mol: input molecule
            mol_center: coordinate of the center of the molecule
            highlight_atom: list of atom indices to be highlighted
            highlight_bond: list of bond indices to be highlighted
        """
        # we offset the drawer such that the image will be drawn at the center
        width, height = self._size_info.get(dm.unique_id(mol))
        self.canvas.SetOffset(int(mol_center.x - width / 2), int(mol_center.y - height / 2))
        self.canvas.SetColour((0, 0, 0, 1))
        self.canvas.DrawMolecule(
            mol,
            legend=self._legends.get(dm.unique_id(mol), ""),
            highlightAtoms=highlight_atom,
            highlightBonds=highlight_bond,
        )
