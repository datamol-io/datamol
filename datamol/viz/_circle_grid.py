from typing import List
from typing import Tuple
from typing import Any
from typing import Optional

import io
import math
import random

from matplotlib.font_manager import FontManager

from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image, ImageDraw, ImageFont, ImageOps


def circle_grid(
    center_mol: Chem.rdchem.Mol,
    circle_mols: List[List[Chem.rdchem.Mol]],
    legend: Optional[str] = None,
    mol_size: Tuple[int, int] = (200, 200),
    circle_margin: int = 50,
    act_mapper: Optional[dict] = None,
):
    """Show molecules in concentric rings, with one molecule at the center

    Args:
        center_mol (Chem.Mol): Molecule at center
        circle_mols (list of list of <Chem.Mol>): List of molecule for each concentric circle around the center mol
        mol_size (tuple, optional): Tuple of width and height for each molecule
        circle_margin (int, optional): Margin between the circle layers
        act_mapper (dict): Map each molecule to a dictionary of activity
    """
    return MolsCircleGrid(center_mol, circle_mols, legend, mol_size, circle_margin, act_mapper)


class MolsCircleGrid:
    def __init__(
        self,
        center_mol: Chem.rdchem.Mol,
        circle_mols: List[List[Chem.rdchem.Mol]],
        legend: Optional[str] = None,
        mol_size: Tuple[int, int] = (200, 200),
        circle_margin: int = 50,
        act_mapper: Optional[dict] = None,
    ):
        """Show molecules in concentric rings, with one molecule at the center

        Args:
            center_mol: Molecule at center
            circle_mols: List of molecule for each concentric circle around the center mol
            mol_size: Tuple of width and height for each molecule
            circle_margin: Margin between the circle layers
            act_mapper: Map each molecule to a dictionary of activity
        """
        self.circle_mols = circle_mols
        self.circle_count = len(self.circle_mols)
        self.legend = legend or ""
        self.margin = circle_margin
        self.center_mol = center_mol
        self.mol_size = mol_size
        size = (max(mol_size) + self.margin) * (self.circle_count + 1)
        self.size = size
        self.image = Image.new(mode="RGBA", size=(size, size), color=(255, 255, 255, 0))
        self.midpoint = size // 2
        self.draw = None
        self.act_mapper = act_mapper or {}
        self._draw()

    def show(self, crop=False):
        if crop:
            crop_img = ImageOps.crop(self.image, border=1)
        else:
            crop_img = self.image
        return crop_img.show()

    def save(self, filename):
        self.image.save(filename)

    def _draw(self):
        """Create circles and slices in-memory"""
        draw = ImageDraw.Draw(self.image)
        self.draw = draw
        all_radius = self._draw_circles(draw)
        self._draw_center_mol()
        self._draw_ring_mols(all_radius)
        font = None
        w, h = draw.textsize(self.legend)
        try:
            fn = FontManager()
            fontpath = fn.findfont("Droid sans")
            font = ImageFont.truetype(fontpath, 12 * self.size // 800)
            w, h = font.getsize(self.legend)
        except:
            pass
        draw.text(
            ((self.size // 2 - w) - 2, self.size - 2 * h),
            self.legend,
            fill="black",
            font=font,
        )
        del draw
        self.draw = None

    def _repr_png_(self):
        bio = io.BytesIO()
        self.image.save(bio, format="PNG")
        return bio.getvalue()

    def _draw_circles(self, draw):
        if self.circle_count <= 0:
            return []
        radius_step = int(self.midpoint / (self.circle_count + 1))
        radius_list = []
        full_range = range(0, self.size // 2, radius_step)

        for i, radius in enumerate(full_range):
            radius += self.margin // 2
            bounding_box = [
                (self.midpoint - radius, self.midpoint - radius),
                (self.midpoint + radius, self.midpoint + radius),
            ]
            if radius > self.margin:
                transp = int(255 - (200 * (i - 1) / len(full_range)))
                draw.arc(bounding_box, 0, 360, fill=(190, 190, 190, transp))
            radius_list.append(radius + radius_step)
        return radius_list

    def _draw_mol_at(
        self,
        mol,
        center_x,
        center_y,
        mol_size=None,
        act_dict={},
        center=False,
        **kwargs: Any,
    ):
        img = mol
        if mol_size is None:
            mol_size = self.mol_size

        if isinstance(mol, Chem.Mol):
            img = Draw.MolToImage(mol, mol_size, kekulize=True, fitImage=True, **kwargs)

        width, height = img.size
        self.image.paste(img, (int(center_x - width / 2), int(center_y - height / 2)))
        txt = []
        for prop, propval in act_dict.items():
            if not isinstance(propval, str):
                propval = "{:.2f}".format(propval)
            txt.append(f"{prop}: {propval}")
        if txt and self.draw is not None:
            txt = "\n".join(txt)
            font = None
            w, h = self.draw.multiline_textsize(txt)
            try:
                fn = FontManager()
                fontpath = fn.findfont("Droid sans")
                font = ImageFont.truetype(fontpath, 18 + center * 8)
                w, h = self.draw.multiline_textsize(txt, font=font)
            except:
                pass

    def _draw_center_mol(self):
        self._draw_mol_at(
            self.center_mol,
            self.midpoint,
            self.midpoint,
            mol_size=[x + self.margin for x in self.mol_size],
            act_dict=self.act_mapper.get(self.center_mol, {}),
            center=True,
        )

    def _draw_ring_mols(self, radius_list):
        for i, mols in enumerate(self.circle_mols):
            radius = radius_list[i]
            ni = len(mols)
            rand_unit = random.random() * 2 * math.pi
            for k, mol in enumerate(mols):
                center_x = radius * math.cos(2 * k * math.pi / ni + rand_unit) + self.midpoint
                center_y = radius * math.sin(2 * k * math.pi / ni + rand_unit) + self.midpoint
                self._draw_mol_at(mol, center_x, center_y, act_dict=self.act_mapper.get(mol, {}))
