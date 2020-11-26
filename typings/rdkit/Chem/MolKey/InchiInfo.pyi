from rdkit import Chem as Chem
from rdkit.Chem import inchi as inchi
from typing import Any

console: Any
UPD_APP: Any
version_re: Any
reconnected_re: Any
fixed_h_re: Any
isotope_re: Any
stereo_re: Any
stereo_all_re: Any
undef_stereo_re: Any
all_stereo_re: Any
defined_stereo_re: Any
h_layer_re: Any
mobile_h_group_re: Any
mobile_h_atoms_re: Any

class InchiInfo:
    parsed_inchi: Any = ...
    def __init__(self, inchi_str: Any) -> None: ...
    def get_sp3_stereo(self): ...
    def get_mobile_h(self): ...
