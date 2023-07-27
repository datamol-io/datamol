# This file is thanks Christian W. Feldman
# Christian Feldman (2021) lassohighlight [sourcecode]. https://github.com/c-feldmann/lassohighlight.

# features to add to this project
# - flag to determine color palette
# - possibility to do this for multiple target molecules at once
# - have the option to write to a file like to_image

from typing import List, Iterator, Tuple, Union, Optional, Any, cast

from collections import defaultdict
from collections import namedtuple

from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.rdmolops import Get3DDistanceMatrix
from rdkit.Geometry.rdGeometry import Point2D

from loguru import logger

import numpy as np
import datamol as dm

from datamol.types import RDKitColor
from datamol.types import DatamolColor

from .utils import drawer_to_image
from .utils import prepare_mol_for_drawing
from .utils import to_rdkit_color


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


def _draw_substructurematch(
    canvas: rdMolDraw2D.MolDraw2D,
    mol: dm.Mol,
    indices: Union[list, str],
    rel_radius: float = 0.3,
    rel_width: float = 0.5,
    line_width: int = 2,
    color: Optional[RDKitColor] = None,
    offset: Optional[Tuple[int, int]] = None,
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
        offset: The offset in raw coordinates for drawing the highlighting given the atom position in the grid.
    """

    prior_lw = canvas.LineWidth()
    canvas.SetLineWidth(line_width)
    canvas.SetFillPolys(False)
    # Setting color
    # #  Default color is gray.
    if not color:
        color = (0.5, 0.5, 0.5, 1)
    canvas.SetColour(tuple(color))

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
    if offset is None:
        offset = Point2D(0, 0)

    for _, at_manager in a_obj_dict.items():
        # A circle is drawn to atoms without outgoing connections
        if not at_manager.bonds:
            pos_list1 = _arch_points(r, 0, np.pi * 2, 60)
            pos_list1[:, 0] += at_manager.pos[0]
            pos_list1[:, 1] += at_manager.pos[1]
            points = [Point2D(*c) for c in pos_list1]
            points = [canvas.GetDrawCoords(p) + offset for p in points]
            canvas.DrawPolygon(points, rawCoords=True)

        # A arch is drawn between attachment points of neighbouring bonds
        for points in at_manager.get_arch_attachment_points():
            pos_list1 = _arch_points(r, points[0], points[1], 20)
            # Translating arch from origin to atom position
            pos_list1[:, 0] += at_manager.pos[0]
            pos_list1[:, 1] += at_manager.pos[1]
            # Transforming points to RDKit Objects
            points = [Point2D(*c) for c in pos_list1]
            points = [canvas.GetDrawCoords(p) + offset for p in points]
            canvas.DrawPolygon(points, rawCoords=True)

        # Drawing lines parallel to each bond
        # EN: we need to convert the coordinates to raw drawing coordinates
        # since RDKit does not expose getDrawTransformers or getAtomCoords
        # and we can't operate in molecule coordinates.
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
            atom_i_left = Point2D(*atom_i_left_at)
            atom_i_left = canvas.GetDrawCoords(atom_i_left) + offset
            atom_j_right = Point2D(*atom_j_right_at)
            atom_j_right = canvas.GetDrawCoords(atom_j_right) + offset

            atom_i_right = Point2D(*atom_i_right_at)
            atom_i_right = canvas.GetDrawCoords(atom_i_right) + offset
            atom_j_left = Point2D(*atom_j_left_at)
            atom_j_left = canvas.GetDrawCoords(atom_j_left) + offset

            canvas.DrawLine(atom_i_left, atom_j_right, rawCoords=True)
            canvas.DrawLine(atom_i_right, atom_j_left, rawCoords=True)
    # restoring prior line width
    canvas.SetLineWidth(prior_lw)


def _draw_multi_matches(
    canvas: rdMolDraw2D.MolDraw2D,
    mol: dm.Mol,
    indices_set_lists: List[Union[list, str]],
    r_min: float = 0.3,
    r_dist: float = 0.13,
    relative_bond_width: float = 0.5,
    color_list: Optional[List[DatamolColor]] = None,
    line_width: int = 2,
    offset: Optional[Tuple[int, int]] = None,
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
        offset: The offset in raw coordinates for drawing the highlighting given the atom position in the grid.
    """
    # If no colors are given, all substructures are depicted in gray.
    if color_list is None or len(color_list) == 0:
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
            color=to_rdkit_color(color),
            line_width=line_width,
            offset=offset,
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
    target_molecules: Union[str, dm.Mol, List[Union[str, dm.Mol]]],
    search_molecules: Union[str, List[str], dm.Mol, List[dm.Mol]] = None,
    atom_indices: Optional[Union[List[int], List[List[int]]]] = None,
    legends: Union[List[Union[str, None]], str, None] = None,
    n_cols: int = 4,
    mol_size: Tuple[int, int] = (300, 300),
    use_svg: Optional[bool] = True,
    draw_mols_same_scale: bool = True,
    r_min: float = 0.3,
    r_dist: float = 0.13,
    relative_bond_width: float = 0.5,
    color_list: Optional[List[DatamolColor]] = None,
    line_width: int = 2,
    scale_padding: float = 1.0,
    verbose: bool = False,
    **kwargs: Any,
):
    """Create an image of a list of molecules with substructure matches using lasso-based highlighting.
    Substructure matching is optional and it's also possible to pass a list of list of atom indices to highlight.

    Args:
        target_molecules:  One or a list of molecules to be highlighted.
        search_molecules: The substructure to be highlighted.
        atom_indices: Atom indices to be highlighted substructure.
        legends: A string or a list of string as legend for every molecules.
        n_cols: Number of molecules per column.
        mol_size: The size of the image to be returned
        use_svg: Whether to return an svg or png image
        draw_mols_same_scale: Whether to draw all the molecules on the same scale. This has the same effect as the `drawMolsSameScale` of the drawing options (which cannot be applied).
        r_min: Radius of the smallest circle around atoms. Length is relative to average bond length (1 = avg. bond len).
        r_dist: Incremental increase of radius for the next substructure.
        relative_bond_width: Distance of line to "bond" (line segment between the two atoms). Size is relative to `atom_radius`.
        color_list: List of tuples with RGBA or RGB values specifying the color of the highlighting.
        line_width: width of drawn lines.
        scale_padding: Padding around the molecule when drawing to scale.
        verbose: Whether to print the verbose information.
        **kwargs: Additional arguments to pass to the drawing function. See RDKit
            documentation related to `MolDrawOptions` for more details at
            https://www.rdkit.org/docs/source/rdkit.Chem.Draw.rdMolDraw2D.html.
    """

    ## Step 0: Input validation

    if search_molecules is None:
        search_molecules = []

    # Make `search_molecules` a list if it is not already
    if not isinstance(search_molecules, (list, tuple)):
        search_molecules = [search_molecules]

    # Convert search molecules into RDKit patterns
    search_molecules = list(search_molecules)
    for i, search_mol in enumerate(search_molecules):
        if isinstance(search_mol, str):
            search_molecules[i] = dm.from_smarts(search_mol)

        if search_molecules[i] is None or not isinstance(search_molecules[i], dm.Mol):
            raise ValueError(
                f"Please enter valid search molecules or smarts: {search_molecules[i]}"
            )

    if not isinstance(target_molecules, (list, tuple)):
        target_molecules = [target_molecules]

    if n_cols is None:
        n_cols = 4

    n_cols = min(n_cols, len(target_molecules))
    n_rows = len(target_molecules) // n_cols
    if len(target_molecules) % n_cols:
        n_rows += 1

    if legends is None:
        legends = [""] * len(target_molecules)
    elif isinstance(legends, str):
        legends = [legends] * len(target_molecules)

    ## Step 1: setup drawer and canvas
    if use_svg:
        drawer = rdMolDraw2D.MolDraw2DSVG(
            mol_size[0] * n_cols,
            mol_size[1] * n_rows,
            mol_size[0],
            mol_size[1],
        )
    else:
        drawer = rdMolDraw2D.MolDraw2DCairo(
            mol_size[0] * n_cols,
            mol_size[1] * n_rows,
            mol_size[0],
            mol_size[1],
        )

    # Setting the drawing options
    draw_options = drawer.drawOptions()
    kwargs = kwargs.copy()
    del_attr = []
    for k, v in kwargs.items():
        if hasattr(draw_options, k):
            setattr(draw_options, k, v)
            del_attr.append(k)
    kwargs = {k: v for k, v in kwargs.items() if k not in del_attr}
    # we will pass the remaining kwargs to the draw function later on

    if color_list is None:
        color_list = DEFAULT_LASSO_COLORS

    mols_to_draw = []
    atoms_idx_list = []
    # scaling parameters
    min_scale_val = Point2D(np.inf, np.inf)
    max_scale_val = Point2D(-np.inf, -np.inf)
    color_list_indexes = []
    ## Step 2: Prepare all the inputs for drawing
    for mol_idx, target_molecule in enumerate(target_molecules):
        # check if the input is valid
        if target_molecule is None:
            raise ValueError("Please enter a valid target molecule or smiles")

        if isinstance(target_molecule, str):
            target_molecule = dm.to_mol(target_molecule)

        # Always make the type checker happy
        target_mol = cast(dm.Mol, target_molecule)
        # Match the search molecules or SMARTS to the target molecule
        atom_idx_list = []
        target_color_list = []
        for color_idx, search_mol in enumerate(search_molecules):
            matches = target_mol.GetSubstructMatches(search_mol)
            if matches:
                matched_atoms = set.union(*[set(x) for x in matches])
                atom_idx_list.append(matched_atoms)
                target_color_list.append(color_idx % len(color_list))
            elif verbose:
                logger.warning(f"No matching substructures found for {dm.to_smarts(search_mol)}")

        ## Add the atom indices to the list if any
        if atom_indices is not None:
            if not isinstance(atom_indices[0], (list, tuple)):
                atom_indices_list_of_list = [atom_indices]
            else:
                atom_indices_list_of_list = atom_indices
            atom_idx_list += atom_indices_list_of_list

        atoms_idx_list.append(atom_idx_list)
        color_list_indexes.append(target_color_list)

        ## Prepare the molecule for drawing and draw it
        mol = prepare_mol_for_drawing(target_molecule, kekulize=True)
        if mol is None:
            raise ValueError(
                f"A molecule {mol_idx} has failed to be prepared by `prepare_mol_for_drawing`."
            )

        conf = mol.GetConformer()
        coordinates = conf.GetPositions()
        min_scale_val.x = min(min_scale_val.x, min(coordinates[:, 0]))
        min_scale_val.y = min(min_scale_val.y, min(coordinates[:, 1]))
        max_scale_val.x = max(max_scale_val.x, max(coordinates[:, 0]))
        max_scale_val.y = max(max_scale_val.y, max(coordinates[:, 1]))
        mols_to_draw.append(mol)

    # Setting up the coordinate system by drawing the molecules as a grid
    # EN: the following is edge-case free after trying 6 different logics, but may break if RDKit changes the way it draws molecules
    scaling_val = Point2D(scale_padding, scale_padding)

    try:
        drawer.DrawMolecules(mols_to_draw, legends=legends, **kwargs)
    except Exception:
        raise ValueError(
            "Failed to draw molecules. Some arguments neither match expected MolDrawOptions, nor DrawMolecule inputs. Please check the input arguments."
        )
    drawer.ClearDrawing()
    if draw_mols_same_scale:
        drawer.SetScale(
            mol_size[0], mol_size[1], min_scale_val - scaling_val, max_scale_val + scaling_val
        )

    for ind, (mol, atom_idx_list) in enumerate(zip(mols_to_draw, atoms_idx_list)):
        h_pos, w_pos = np.unravel_index(ind, (n_rows, n_cols))
        offset_x = int(w_pos * mol_size[0])
        offset_y = int(h_pos * mol_size[1])
        drawer.SetOffset(offset_x, offset_y)
        drawer.DrawMolecule(mol, legend=legends[ind], **kwargs)
        offset = None
        if draw_mols_same_scale:
            offset = drawer.Offset()
            # EN: if the molecule has a legend we need to offset the highlight by the height of the legend
            # we also need to account for the padding around the molecule
            if legends[ind]:
                # geometry is hard
                padding_fraction = (drawer.drawOptions().legendFraction * mol_size[1]) // 2
                offset -= Point2D(0, padding_fraction)

        if len(atom_idx_list) > 0:
            dm.viz._lasso_highlight._draw_multi_matches(
                drawer,
                mol,
                atom_idx_list,
                r_min=r_min,
                r_dist=r_dist,
                relative_bond_width=relative_bond_width,
                line_width=line_width,
                color_list=[color_list[i] for i in color_list_indexes[ind]],
                offset=offset,
            )
        elif verbose:
            logger.warning("No matches found for the given search molecules")

    drawer.FinishDrawing()
    # NOTE(hadim): process the drawer object to return the image type matching the same behavior as RDkit and `datamol.to_image()`
    return drawer_to_image(drawer)
