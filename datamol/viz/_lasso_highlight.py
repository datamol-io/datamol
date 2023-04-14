# This file is thanks Christian W. Feldman
# Christian Feldman (2021) lassohighlight [sourcecode]. https://github.com/c-feldmann/lassohighlight.

# features to add to this project
# - flag to determine color palette
# - possibility to do this for multiple target molecules at once
# - have the option to write to a file like to_image

from typing import List, Iterator, Tuple, Union, Optional, Any

from collections import defaultdict
from collections import namedtuple

from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.rdmolops import Get3DDistanceMatrix
from rdkit.Geometry.rdGeometry import Point2D

from loguru import logger

import numpy as np
import datamol as dm

from .utils import drawer_to_image
from .utils import prepare_mol_for_drawing


def _angle_to_coord(center: np.ndarray, angle: float, radius: float) -> np.ndarray:
    """Determines a point relative to the center with distance (radius) at given angle.
    Angles are given in rad and 0 rad correspond to north of the center point.

    args:
        center: The center point.
        angle: The angle in rad.
        radius: The distance of the point to the center.
    """

    x = radius * np.sin(angle)
    y = radius * np.cos(angle)
    x += center[0]
    y += center[1]
    return np.array([x, y])


def _arch_points(radius: float, start_ang: float, end_ang: float, n: int) -> np.ndarray:
    """Returns an array of the shape (2, n) with equidistant points on the arch defined by
    given radius and angles. Angles are given in rad.

    args:
        radius: The radius of the arch.
        start_ang: The start angle of the arch.
        end_ang: The end angle of the arch.
        n: The number of points to return.
    """
    angles = np.linspace(start_ang, end_ang, n)
    x = radius * np.sin(angles)
    y = radius * np.cos(angles)
    return np.vstack([x, y]).T


def _angle_between(center: np.ndarray, pos: np.ndarray) -> np.ndarray:
    """Calculates the angle in rad between two points.
    An angle of 0 corresponds to north of the center.

    args:
        center: The center point.
        pos: The point to calculate the angle to.
    """

    diff = pos - center
    return np.arctan2(diff[0], diff[1])


def _avg_bondlen(mol: dm.Mol) -> float:
    """Calculates the average bond length of an dm.Mol object.

    args:
        mol: The dm.Mol object.
    """
    distance_matrix = Get3DDistanceMatrix(mol)

    bondlength_list = []
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()
        bondlength_list.append(distance_matrix[a1, a2])

    return float(np.mean(bondlength_list))


Bond = namedtuple("Bond", ["angle", "neighbour_id", "bond_id"])


class _AttachmentPointManager:
    """AnchorManager is an invisible overlay for RDKit Atoms storing positions for
    arches and bond-attachment-points.
    """

    def __init__(self, position: np.ndarray, radius: float, bond_width: float):
        self.pos = position
        self.bond_width = bond_width
        self.radius = radius
        self.bonds: List[Bond] = []
        self.bond_attachment_points: Optional[dict] = None

    def add_bond(self, angle: float, neighbor_id: int, bond_id: int):
        self.bonds.append(Bond(angle, neighbor_id, bond_id))

    @property
    def delta_angle(self) -> np.ndarray:
        return np.arcsin(self.bond_width / self.radius)

    def generate_attachment_points(self):
        sorted_bonds = sorted(self.bonds, key=lambda x: x.angle)
        self.bond_attachment_points = dict()
        for i, bond in enumerate(sorted_bonds):
            alpha_left = bond.angle - self.delta_angle
            alpha_right = bond.angle + self.delta_angle

            d_left = self.radius
            d_right = self.radius

            # Handling intersecting bonds
            # # Intersections with previous bonds
            if i == 0:  # For first bond the previous bond is the last bond. Subtracting 2 pi.
                prev_bond_angle = sorted_bonds[-1].angle - np.pi * 2
            else:
                prev_bond_angle = sorted_bonds[i - 1].angle

            # #  If both points intersect the mean angle is calculated.
            if prev_bond_angle + self.delta_angle >= alpha_left:
                alpha_left = np.mean([prev_bond_angle + self.delta_angle, alpha_left])

                a_rhombus = bond.angle - prev_bond_angle

                len_rhombus = self.bond_width / np.sin(a_rhombus)
                # Radius is altered to match the intersecting position
                d_left = 2 * len_rhombus * np.cos(a_rhombus / 2)

            # # Intersections with following bonds
            if i + 1 == len(sorted_bonds):
                next_bond_angle = sorted_bonds[0].angle + np.pi * 2
            else:
                next_bond_angle = sorted_bonds[i + 1].angle

            if next_bond_angle - self.delta_angle <= alpha_right:
                alpha_right = np.mean([next_bond_angle - self.delta_angle, alpha_right])

                a_rhombus = next_bond_angle - bond.angle
                len_rhombus = self.bond_width / np.sin(a_rhombus)
                d_right = 2 * len_rhombus * np.cos(a_rhombus / 2)

            self.bond_attachment_points[bond.bond_id] = [
                (alpha_left, d_left),
                (alpha_right, d_right),
            ]
        return self

    def get_arch_attachment_points(self) -> Iterator[Tuple[float, float]]:
        """Points between bonds which are drawn as arch."""

        if self.bond_attachment_points is None:
            raise ValueError(
                "Attachment points have to be generated first with `generate_attachment_points()`"
            )

        if self.bonds:
            sorted_bonds = sorted(self.bonds, key=lambda x: x[0])
            _, _, bond_keys = zip(*sorted_bonds)
            for i, k in enumerate(bond_keys):
                if i == 0:
                    start_angle = self.bond_attachment_points[bond_keys[-1]][1][0] - np.pi * 2
                else:
                    start_angle = self.bond_attachment_points[bond_keys[i - 1]][1][0]
                end_angle = self.bond_attachment_points[k][0][0]
                if np.isclose(start_angle % (np.pi * 2), end_angle % (np.pi * 2)):
                    continue
                yield start_angle, end_angle


ColorTuple = Union[Tuple[float, float, float, float], Tuple[float, float, float]]


def _draw_substructurematch(
    canvas: rdMolDraw2D.MolDraw2D,
    mol: dm.Mol,
    indices: Union[list, str],
    rel_radius: float = 0.3,
    rel_width: float = 0.5,
    line_width: int = 2,
    color: Optional[ColorTuple] = None,
) -> None:
    """Draws the substructure defined by (atom-) `indices`, as lasso-highlight onto `canvas`.

    args:
        canvas : RDKit Canvas, where highlighting is drawn to.
        mol: Atoms from the molecule `mol` are takes as positional reference for the highlighting.
        indices: Atom indices delineating highlighted substructure.
        rel_radius: Radius of the circle around atoms. Length is relative to average bond length (1 = avg. bond len).
        rel_width: Distance of line to "bond" (line segment between the two atoms). Size is relative to `atom_radius`.
        line_width: width of drawn lines.
        color: Tuple with RGBA or RGB values specifying the color of the highlighting.
    """

    prior_lw = canvas.LineWidth()
    canvas.SetLineWidth(line_width)
    canvas.SetFillPolys(False)
    # Setting color
    # #  Default color is gray.
    if not color:
        color = (0.5, 0.5, 0.5, 1)
    canvas.SetColour(color)

    # Selects first conformer and calculates the mean bond length
    conf = mol.GetConformer(0)
    avg_len = _avg_bondlen(mol)
    r = avg_len * rel_radius
    w = r * rel_width

    a_obj_dict = dict()  # Dictionary for atoms delineating highlighted substructure.
    for atom in mol.GetAtoms():
        a_idx = atom.GetIdx()
        if a_idx not in indices:
            continue

        # 2D-coordinates of Atom
        atom_pos = conf.GetAtomPosition(a_idx)
        atom_pos = np.array([atom_pos.x, atom_pos.y])

        # Initializing an AttachmentPointManager centered at the atom position
        at_manager = _AttachmentPointManager(atom_pos, r, w)

        # Adding Bonds to the AttachmentPointManager
        for bond in atom.GetBonds():
            bond_atom1 = bond.GetBeginAtomIdx()
            bond_atom2 = bond.GetEndAtomIdx()
            neigbor_idx = bond_atom1 if bond_atom2 == a_idx else bond_atom2
            if neigbor_idx not in indices:
                continue
            neigbor_pos = conf.GetAtomPosition(neigbor_idx)
            neigbor_pos = np.array([neigbor_pos.x, neigbor_pos.y])
            bond_angle = float(_angle_between(atom_pos, neigbor_pos))
            bond_angle = bond_angle % (2 * np.pi)  # Assuring 0 <= bond_angle <= 2 pi
            at_manager.add_bond(bond_angle, neigbor_idx, bond.GetIdx())
        at_manager.generate_attachment_points()
        a_obj_dict[a_idx] = at_manager

    added_bonds = set()
    for idx, at_manager in a_obj_dict.items():
        # A circle is drawn to atoms without outgoing connections
        if not at_manager.bonds:
            pos_list1 = _arch_points(r, 0, np.pi * 2, 60)
            pos_list1[:, 0] += at_manager.pos[0]
            pos_list1[:, 1] += at_manager.pos[1]
            points = [Point2D(*c) for c in pos_list1]
            canvas.DrawPolygon(points)

        # A arch is drawn between attachment points of neighbouring bonds
        for points in at_manager.get_arch_attachment_points():
            pos_list1 = _arch_points(r, points[0], points[1], 20)
            # Translating arch from origin to atom position
            pos_list1[:, 0] += at_manager.pos[0]
            pos_list1[:, 1] += at_manager.pos[1]
            # Transforming points to RDKit Objects
            points = [Point2D(*c) for c in pos_list1]
            canvas.DrawPolygon(points)

        # Drawing lines parallel to each bond
        for bond in at_manager.bonds:
            if bond.bond_id in added_bonds:
                continue
            added_bonds.add(bond.bond_id)
            bnd_points = at_manager.bond_attachment_points[bond.bond_id]

            atom_i_left_at = _angle_to_coord(at_manager.pos, *bnd_points[0])
            atom_i_right_at = _angle_to_coord(at_manager.pos, *bnd_points[1])
            atom_j = a_obj_dict[bond.neighbour_id]
            atom_j_left_at = _angle_to_coord(
                atom_j.pos, *atom_j.bond_attachment_points[bond.bond_id][0]
            )
            atom_j_right_at = _angle_to_coord(
                atom_j.pos, *atom_j.bond_attachment_points[bond.bond_id][1]
            )
            canvas.DrawLine(Point2D(*atom_i_left_at), Point2D(*atom_j_right_at))
            canvas.DrawLine(Point2D(*atom_i_right_at), Point2D(*atom_j_left_at))
    # restoring prior line width
    canvas.SetLineWidth(prior_lw)


def _draw_multi_matches(
    canvas: rdMolDraw2D.MolDraw2D,
    mol: dm.Mol,
    indices_set_lists: List[Union[list, str]],
    r_min: float = 0.3,
    r_dist: float = 0.13,
    relative_bond_width: float = 0.5,
    color_list: Optional[List[ColorTuple]] = None,
    line_width: int = 2,
):
    """Draws multiple substructure matches on a canvas.

    args:
        canvas : RDKit Canvas, where highlighting is drawn to.
        mol: Atoms from the molecule `mol` are takes as positional reference for the highlighting.
        indices_set_lists: Atom indices delineating highlighted substructure.
        r_min: Radius of the smallest circle around atoms. Length is relative to average bond length (1 = avg. bond len).
        r_dist: Incremental increase of radius for the next substructure.
        relative_bond_width: Distance of line to "bond" (line segment between the two atoms). Size is relative to `atom_radius`.
        line_width: width of drawn lines.
        color_list: List of tuples with RGBA or RGB values specifying the color of the highlighting.
    """
    # If no colors are given, all substructures are depicted in gray.
    if color_list is None:
        _color_list = [(0.5, 0.5, 0.5)] * len(indices_set_lists)
    else:
        _color_list = color_list

    if len(_color_list) < len(indices_set_lists):
        colors_to_add = []
        for i in range(len(indices_set_lists) - len(_color_list)):
            colors_to_add.append(_color_list[i % len(_color_list)])
        _color_list.extend(colors_to_add)

    level_manager = defaultdict(set)
    for match_atoms, color in zip(indices_set_lists, _color_list):
        used_levels = set.union(*[level_manager[a] for a in match_atoms])

        if len(used_levels) == 0:
            free_levels = {0}
        else:
            max_level = max(used_levels)
            free_levels = set(range(max_level)) - used_levels

        if free_levels:
            draw_level = min(free_levels)
        else:
            draw_level = max(used_levels) + 1

        for a in match_atoms:
            level_manager[a].add(draw_level)

        ar = r_min + r_dist * draw_level
        _draw_substructurematch(
            canvas,
            mol,
            match_atoms,
            rel_radius=ar,
            rel_width=max(relative_bond_width, ar),
            color=color,
            line_width=line_width,
        )


DEFAULT_LASSO_COLORS = [
    (1, 0, 0, 1),  # red
    (0, 0.5, 1, 1),  # blue
    (1, 0.5, 0, 1),  # orange
    (0, 1, 0, 1),  # green
    (1, 1, 0, 1),  # yellow
    (0, 0, 0.5, 1),  # dark blue
]


def lasso_highlight_image(
    target_molecule: Union[str, dm.Mol],
    search_molecules: Union[str, List[str], dm.Mol, List[dm.Mol]],
    mol_size: Tuple[int, int] = (300, 300),
    use_svg: Optional[bool] = True,
    r_min: float = 0.3,
    r_dist: float = 0.13,
    relative_bond_width: float = 0.5,
    color_list: Optional[List[ColorTuple]] = None,
    line_width: int = 2,
    **kwargs: Any,
):
    """Create an image of a molecule with substructure matches using lasso-based highlighting.

    args:
        target_molecule: The molecule to be highlighted
        search_molecules: The substructure to be identified
        mol_size: The size of the image to be returned
        use_svg: Whether to return an svg or png image
        r_min: Radius of the smallest circle around atoms. Length is relative to average bond length (1 = avg. bond len).
        r_dist: Incremental increase of radius for the next substructure.
        relative_bond_width: Distance of line to "bond" (line segment between the two atoms). Size is relative to `atom_radius`.
        line_width: width of drawn lines.
        color_list: List of tuples with RGBA or RGB values specifying the color of the highlighting.
        **kwargs: Additional arguments to pass to the drawing function. See RDKit
            documentation related to `MolDrawOptions` for more details at
            https://www.rdkit.org/docs/source/rdkit.Chem.Draw.rdMolDraw2D.html.
    """

    # check if the input is valid
    if target_molecule is None or (isinstance(target_molecule, str) and len(target_molecule) == 0):
        raise ValueError("Please enter a valid target molecule or smiles")

    if search_molecules is None or (
        isinstance(search_molecules, str) and len(search_molecules) == 0
    ):
        raise ValueError("Please enter valid search molecules or smarts")

    # less than 1 throws File parsing error: PNG header not recognized over 5,000 leads to a DecompressionBombError later on
    if mol_size[0] < 1 or mol_size[0] > 5000 or mol_size[1] < 1 or mol_size[1] > 5000:
        raise ValueError(
            "To avoid errors please choose a number between 1-5000 for the canvas width or height"
        )

    if isinstance(target_molecule, str):
        target_molecule = dm.to_mol(target_molecule)

    mol = prepare_mol_for_drawing(target_molecule, kekulize=True)

    if mol is None:
        raise ValueError("The molecule has failed to be prepared by `prepare_mol_for_drawing`.")

    if use_svg:
        drawer = rdMolDraw2D.MolDraw2DSVG(mol_size[0], mol_size[1])
    else:
        drawer = rdMolDraw2D.MolDraw2DCairo(mol_size[0], mol_size[1])

    # Setting the drawing options
    draw_options = drawer.drawOptions()
    for k, v in kwargs.items():
        if not hasattr(draw_options, k):
            raise ValueError(
                f"Invalid drawing option: {k}={v}. Check `rdkit.Chem.Draw.rdMolDraw2D.MolDrawOptions` for valid ones."
            )
        else:
            setattr(draw_options, k, v)

    # Setting up the coordinate system by drawing and erasing molecule
    drawer.DrawMolecule(mol)
    drawer.ClearDrawing()

    # get the atom indices for the search molecules
    atom_idx_list = []
    if isinstance(search_molecules, str):
        smart_obj = dm.to_mol(search_molecules)
        matches = mol.GetSubstructMatches(smart_obj)
        if not matches:
            logger.warning(f"no matching substructure found for {search_molecules}")
        else:
            matched_atoms = set.union(*[set(x) for x in matches])
            atom_idx_list.append(matched_atoms)

    elif isinstance(search_molecules, dm.Mol):
        matches = mol.GetSubstructMatches(search_molecules)
        if not matches:
            logger.warning(f"no matching substructure found for {dm.to_smiles(search_molecules)}")
        else:
            matched_atoms = set.union(*[set(x) for x in matches])
            atom_idx_list.append(matched_atoms)

    elif len(search_molecules) and isinstance(search_molecules[0], str):
        for smart_str in search_molecules:
            smart_obj = dm.to_mol(smart_str)
            matches = mol.GetSubstructMatches(smart_obj)
            if not matches:
                logger.warning(f"no matching substructure found for {smart_str}")
            else:
                matched_atoms = set.union(*[set(x) for x in matches])
                atom_idx_list.append(matched_atoms)

    elif len(search_molecules) and isinstance(search_molecules[0], dm.Mol):
        for smart_obj in search_molecules:
            matches = mol.GetSubstructMatches(smart_obj)
            if not matches:
                logger.warning(f"no matching substructure found for {dm.to_smiles(smart_obj)}")
            else:
                matched_atoms = set.union(*[set(x) for x in matches])
                atom_idx_list.append(matched_atoms)

    if color_list is None:
        color_list = DEFAULT_LASSO_COLORS

    if len(atom_idx_list) == 0:
        logger.warning("No matches found for the given search molecules")
    else:
        _draw_multi_matches(
            drawer,
            mol,
            atom_idx_list,
            r_min=r_min,
            r_dist=r_dist,
            relative_bond_width=relative_bond_width,
            line_width=line_width,
            color_list=color_list,
        )

    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()

    # NOTE(hadim): process the drawer object to return the image type matching the same behavior as RDkit and `datamol.to_image()`
    return drawer_to_image(drawer)
