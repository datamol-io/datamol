#This file is thanks Christian W. Feldman
# Christian Feldman (2021) lassohighlight [sourcecode]. https://github.com/c-feldmann/lassohighlight.
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.rdmolops import Get3DDistanceMatrix
from rdkit import Chem
from rdkit.Geometry.rdGeometry import Point2D
import numpy as np
from collections import defaultdict
from collections import namedtuple
from typing import *
from PIL import Image
import io

import datamol as dm
from .utils import prepare_mol_for_drawing

def _angle_to_coord(center: np.ndarray, angle: np.float64, radius: np.float64) -> np.ndarray:
    """Determines a point relative to the center with distance (radius) at given angle.
    Angles are given in rad and 0 rad correspond to north of the center point.
    """
    
    x = radius * np.sin(angle)
    y = radius * np.cos(angle)
    x += center[0]
    y += center[1]
    return np.array([x, y])

def _arch_points(radius: np.float64, start_ang: np.float64, end_ang:np.float64, n: int) -> np.ndarray:
    """Returns an array of the shape (2, n) with equidistant points on the arch defined by 
    given radius and angles. Angles are given in rad.
    """
    angles = np.linspace(start_ang, end_ang, n)
    x = radius * np.sin(angles)
    y = radius * np.cos(angles)
    return np.vstack([x, y]).T

def _angle_between(center: np.ndarray, pos: np.ndarray) -> np.ndarray:
    """Calculates the angle in rad between two points.
    An angle of 0 corresponds to north of the center.
    """
    diff = pos - center
    return np.arctan2(diff[0], diff[1])

def _avg_bondlen(mol: Chem.Mol) -> np.ndarray:
    """Calculates the average bond length of an rdkit.Chem.Mol object.
    """
    distance_matrix = Get3DDistanceMatrix(mol)
    bondlength_list: List = []
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()
        bondlength_list.append(distance_matrix[a1, a2])
    return np.mean(bondlength_list)

Bond = namedtuple("Bond", ["angle", "neighbour_id", "bond_id"])


class _AttachmentPointManager:
    """AnchorManager is an invisible overlay for RDKit Atoms storing positions for arches and bond-attachment-points.
    """

    def __init__(self, position: np.array, radius: np.float64, bond_width: np.float64 ):
        self.pos = position
        self.bond_width = bond_width
        self.radius = radius
        self.bonds: List[Bond] = []
        self.bond_attachment_points: Optional[List[Tuple]] = None

    def add_bond(self, angle: np.float64, neighbor_id: int, bond_id: int):
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

            self.bond_attachment_points[bond.bond_id] = [(alpha_left, d_left), (alpha_right, d_right)]
        return self

    def get_arch_attachment_points(self) -> Iterator[Tuple[float, float]]:
        """Points between bonds which are drawn as arch.
        """
        
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
    mol: Chem.Mol, 
    indices: Union[list, str], 
    rel_radius: float=0.3, 
    rel_width: float=0.5, 
    line_width: int=2, 
    color: ColorTuple=None
    ) -> None :
    """ Draws the substructure defined by (atom-) `indices`, as lasso-highlight onto `canvas`.
    Parameters
    ----------
    canvas : rdMolDraw2D.MolDraw2D
        RDKit Canvas, where highlighting is drawn to.
    mol: Chem.Mol
        Atoms from the molecule `mol` are takes as positional reference for the highlighting.
    indices: Union[list, str]
        Atom indices delineating highlighted substructure.
    rel_radius: float
        Radius of the circle around atoms. Length is relative to average bond length (1 = avg. bond len).
    rel_width: float
        Distance of line to "bond" (line segment between the two atoms). Size is relative to `atom_radius`.
    line_width: int
        width of drawn lines.
    color: ColorTuple
           Tuple with RGBA or RGB values specifying the color of the highlighting.
    Returns
    -------
    None
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
            bond_angle = _angle_between(atom_pos, neigbor_pos)
            bond_angle = bond_angle % (2*np.pi)  # Assuring 0 <= bond_angle <= 2 pi
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
            atom_j_left_at = _angle_to_coord(atom_j.pos, *atom_j.bond_attachment_points[bond.bond_id][0])
            atom_j_right_at = _angle_to_coord(atom_j.pos, *atom_j.bond_attachment_points[bond.bond_id][1])
            canvas.DrawLine(Point2D(*atom_i_left_at), Point2D(*atom_j_right_at))
            canvas.DrawLine(Point2D(*atom_i_right_at), Point2D(*atom_j_left_at))
    # restoring prior line width
    canvas.SetLineWidth(prior_lw)

#TODO switch this over to other doc string
def _draw_multi_matches(
    canvas: rdMolDraw2D.MolDraw2D, 
    mol: Chem.Mol, 
    indices_set_lists: List[Union[list, str]], 
    r_min: float=0.3, 
    r_dist: float=0.13, 
    relative_bond_width: float=0.5, 
    color_list: List[ColorTuple]=None,
    line_width: int=2
    ) -> None:
    """
    Parameters
    ----------
    canvas : rdMolDraw2D.MolDraw2D
        RDKit Canvas, where highlighting is drawn to.
    mol: Chem.Mol
        Atoms from the molecule `mol` are takes as positional reference for the highlighting.
    indices_set_lists: List[Union[list, str]]
        Atom indices delineating highlighted substructure.
    r_min: float
        Radius of the smallest circle around atoms. Length is relative to average bond length (1 = avg. bond len).
    r_dist: float
        Incremental increase of radius for the next substructure.
    relative_bond_width: float
        Distance of line to "bond" (line segment between the two atoms). Size is relative to `atom_radius`.
    line_width: int
        width of drawn lines.
    color_list: List[ColorTuple]
        List of tuples with RGBA or RGB values specifying the color of the highlighting.
    Returns
    -------
    None
    Returns
    -------
    """
    # If no colors are given, all substructures are depicted in gray.
    if color_list is None:
        color_list = [(0.5, 0.5, 0.5)] * len(indices_set_lists)
    if len(color_list) < len(indices_set_lists):
        raise ValueError("Not enough colors for referenced substructures!")

    level_manager = defaultdict(set)
    for match_atoms, color in zip(indices_set_lists, color_list):
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
        _draw_substructurematch(canvas,
                                mol,
                                match_atoms,
                                rel_radius=ar,
                                rel_width=max(relative_bond_width, ar),
                                color=color,
                                line_width=line_width)

def lasso_highlight_image(
    canvas_width: int, 
    canvas_height: int,
    target_molecule: Union[str, dm.Mol],
    search_molecules: Union[str, List[str], dm.Mol, List[dm.Mol]]
    ) -> Image:
    """A generalized interface to access both highlighting options whether the input is as a smiles, smarts or mol

    Args:
        canvas_width (int): _description_
        canvas_height (int): _description_
        target_molecule (Union[str, dm.Mol]): 
        search_molecules (Union[str, List[str], dm.Mol, List[dm.Mol]]): Various ways to enter but the search molecules must be entered as smart
    """
    
    color_dict = { "red":   (1,   0,    0,  1),
            "blue":   (0,   0.5,  1,    1), 
            "orange": (1,   0.5,  0,    1),
            "green": (0,   1, 0, 1),
            "yellow": (1, 1, 0, 1),
            "dark blue": (0, 0, 0.5 , 1)
            }
    
    #less than 1 throws File parsing error: PNG header not recognized over 5,000 leads to a DecompressionBombError later on
    if canvas_width < 1 or canvas_width > 5000 or canvas_height < 1 or canvas_height > 5000:
        raise ValueError(f"To avoid errors please choose a number between 1-5000 for the canvas width or height")
    
    if isinstance(target_molecule, str):
            target_molecule = dm.to_mol(target_molecule)
    
    mol = prepare_mol_for_drawing(target_molecule, kekulize=True)
    d = rdMolDraw2D.MolDraw2DCairo(500, 500)
    d.drawOptions().updateAtomPalette({i: (0, 0, 0, 1) for i in range(100)})
    d.DrawMolecule(mol)
    d.ClearDrawing()

    atom_idx_list = []
        
    if isinstance(search_molecules, str):
        smart_obj = Chem.MolFromSmarts(search_molecules)
        matched_atoms = set.union(*[set(x) for x in mol.GetSubstructMatches(smart_obj)])
        if len(matched_atoms) >=1:
            atom_idx_list.append(matched_atoms)
        
    elif isinstance(search_molecules[0], str):
        for smart_str in search_molecules:
            smart_obj = Chem.MolFromSmarts(smart_str)
            matched_atoms = set.union(*[set(x) for x in mol.GetSubstructMatches(smart_obj)])
            if len(matched_atoms) >=1:
                atom_idx_list.append(matched_atoms)
        
    _draw_multi_matches(d, mol, atom_idx_list, r_min=0.3, r_dist=0.12, relative_bond_width=0.5, color_list=color_dict.values(), line_width=2)
        
    d.DrawMolecule(mol)
    d.FinishDrawing()

    return Image.open(io.BytesIO(d.GetDrawingText()))
        