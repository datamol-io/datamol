from typing import Optional
import collections
from rdkit import Chem
from rdkit.Chem import rdChemReactions
import datamol as dm


IsomerReaction = collections.namedtuple(
    "IsomerReaction",
    ["name", "smarts", "reverse", "acyclic", "triplebond", "doublebond", "use"],
)

ISOMERIZERS = dict(
    R1=IsomerReaction(
        "R1",
        "([*:1]-[*:2].[*:3]-[*:4])>>([*:1]-[*:4].[*:2]-[*:3])",
        "R1",
        True,
        False,
        False,
        True,
    ),
    R1a=IsomerReaction(
        "R1a",
        "([$([*]-[*]):1]-[$([*]-[*:1]):2].[$([*]-[*]):3]-[$([*]-[*:3]):4])>>([*:1]-[*:4].[*:3]-[*:2])",
        None,
        True,
        False,
        False,
        False,
    ),
    R1b=IsomerReaction(
        "R1b",
        "([*:2]-[*:1]-[*:3]-[*:4])>>([*:2]-[*:3]-[*:1]-[*:4])",
        "R1b",
        True,
        False,
        False,
        False,
    ),
    R1c=IsomerReaction(
        "R1c",
        "([*:1]-[*:2]-[*:3]-[*:4])>>([*:4]-[*:2]-[*:3]-[*:1])",
        "R1c",
        True,
        False,
        False,
        False,
    ),
    R2=IsomerReaction(
        "R2",
        "([*:1]=[*:2].[*:3]-[*:4])>>([*:3]-[*:2]-[*:1]-[*:4])",
        "R4",
        False,
        False,
        True,
        True,
    ),
    R3=IsomerReaction(
        "R3",
        "([*:1]=[*:2].[*:3]=[*:4])>>([*:1]-1-[*:2]-[*:3]-[*:4]-1)",
        "R7",
        False,
        False,
        True,
        True,
    ),
    R4=IsomerReaction(
        "R4",
        "([*:2]-[*:1]-[*:4]-[*:3])>>([*:2]-[*:3].[*:1]=[*:4])",
        "R2",
        False,
        False,
        True,
        True,
    ),
    R5=IsomerReaction(
        "R5",
        "([*:3]-[*:4]-[*:1]=[*:2])>>([*:3]-[*:2]-[*:1]=[*:4])",
        "R5",
        True,
        False,
        True,
        True,
    ),
    R6=IsomerReaction(
        "R6",
        "([*:2]=[*:1]-[*:4]=[*:3])>>([*:2]-1-[*:3]-[*:4]=[*:1]-1)",
        "R8",
        False,
        False,
        True,
        True,
    ),
    R7=IsomerReaction(
        "R7",
        "([*:1]-1-[*:2]-[*:3]-[*:4]-1)>>([*:2]=[*:3].[*:1]=[*:4])",
        "R5",
        False,
        False,
        True,
        True,
    ),
    R8=IsomerReaction(
        "R8",
        "([*:3]-1-[*:4]-[*:1]=[*:2]-1)>>([*:3]=[*:2]-[*:1]=[*:4])",
        "R6",
        False,
        False,
        True,
        True,
    ),
    R9=IsomerReaction(
        "R9",
        "([*:1]-1=[*:2]-[*:3]=[*:4]-1)>>([*:1]-1=[*:4]-[*:3]=[*:2]-1)",
        "R9",
        True,
        False,
        True,
        True,
    ),
    R10=IsomerReaction(
        "R10",
        "([*:2]-[*:1]=[*:4]-[*:3])>>([*:2]-[*:3].[*:1]#[*:4])",
        "R15",
        False,
        True,
        True,
        True,
    ),
    R11=IsomerReaction(
        "R11",
        "([*:3]-[*:4]=[*:1]=[*:2])>>([*:3]-[*:2]-[*:1]#[*:4])",
        "R18",
        True,
        True,
        True,
        True,
    ),
    R12=IsomerReaction(
        "R12",
        "([*:2]=[*:1]=[*:4]=[*:3])>>([*:2]-1-[*:3]-[*:4]#[*:1]-1)",
        "R19",
        False,
        True,
        True,
        True,
    ),
    R13=IsomerReaction(
        "R13",
        "([*:2]-1-[*:3]-[*:4]=[*:1]-1)>>([*:2]=[*:3].[*:1]#[*:4])",
        "R16",
        False,
        True,
        True,
        True,
    ),
    R14=IsomerReaction(
        "R14",
        "([*:1]-1=[*:4]-[*:3]=[*:2]-1)>>([*:2]#[*:3].[*:1]#[*:4])",
        "R17",
        False,
        True,
        True,
        True,
    ),
    R15=IsomerReaction(
        "R15",
        "([*:1]#[*:2].[*:3]-[*:4])>>([*:3]-[*:2]=[*:1]-[*:4])",
        "R10",
        False,
        True,
        True,
        True,
    ),
    R16=IsomerReaction(
        "R16",
        "([*:1]#[*:2].[*:3]=[*:4])>>([*:3]-1-[*:4]-[*:1]=[*:2]-1)",
        "R13",
        False,
        True,
        True,
        True,
    ),
    R17=IsomerReaction(
        "R17",
        "([*:1]#[*:2].[*:3]#[*:4])>>([*:1]-1=[*:2]-[*:3]=[*:4]-1)",
        "R14",
        False,
        True,
        True,
        True,
    ),
    R18=IsomerReaction(
        "R18",
        "([*:3]-[*:4]-[*:1]#[*:2])>>([*:3]-[*:2]=[*:1]=[*:4])",
        "R11",
        True,
        True,
        True,
        True,
    ),
    R19=IsomerReaction(
        "R19",
        "([*:3]-1-[*:4]-[*:1]#[*:2]-1)>>([*:3]=[*:2]=[*:1]=[*:4])",
        "R12",
        False,
        True,
        True,
        True,
    ),
)


class IsomerEnumerator:
    """
    Implementation of the isomer enumeration algorithm described in  https://doi.org/10.1186/s13321-017-0252-9

    ..note::
        Due to the chemical reaction used, atom mapping will not be preserved !

    """

    def __init__(
        self,
        allow_cycle: bool = True,
        allow_double_bond: bool = False,
        allow_triple_bond: bool = False,
        rxn_list: Optional[list] = None,
    ):
        """Structural isomer enumeration

        Args:
            allow_cycle (bool, optional): whether to allow transformation involving cycle creation
            allow_double_bond (bool, optional): whether to allow transformation involving at least one double bond
            allow_triple_bond (bool, optional): whether to allow transformation involving at least one triple bond
            rxn_list (list, optional): List of str to specifiy the reactions that should be used
        """
        self.allow_cycle = allow_cycle
        self.allow_triple_bond = allow_triple_bond
        self.allow_double_bond = allow_double_bond
        self.rxn_cache = self._resolve(rxn_list)
        self._ring_smarts = dm.from_smarts("[R]")
        # there should be a lot more of this
        self._impossible_chemistry = [dm.from_smarts("[$([*;r3,r4](=*)=*)]")]
        if not self.rxn_cache:
            raise ValueError("No isomeric transformation matches your")

    def _resolve(self, rxn_list: list):
        """Resolve map of rxn to use

        Args:
            rxn_list: list of reaction to use
        """
        if rxn_list is not None and len(rxn_list) > 0:
            rxn_list = [_rxn for _rxn in rxn_list if _rxn in ISOMERIZERS]
        else:
            rxn_list = list([x for x, val in ISOMERIZERS.items() if val.use])

        final_rxns = dict()

        for x in rxn_list:
            rxn_trans = ISOMERIZERS[x]
            # check for triple bonds
            can_add = (not rxn_trans.triplebond) or self.allow_triple_bond
            # check for double bonds
            can_add &= not (rxn_trans.doublebond) or self.allow_double_bond
            # check for cycles
            can_add &= rxn_trans.acyclic or self.allow_cycle
            if can_add:
                as_rxn = Chem.rdChemReactions.ReactionFromSmarts(rxn_trans.smarts)
                final_rxns[x] = as_rxn
        return final_rxns

    def _clean(self, mol: dm.Mol):
        """Clean and sanitize a molecule during enumeration

        Args:
            mol: molecule to clean
        """
        clean_mol = None
        clean_smiles = ""
        try:
            mol = dm.sanitize_mol(mol)
            Chem.rdmolops.AssignRadicals(mol)
            for a in mol.GetAtoms():
                a.SetNoImplicit(True)
                a.SetNumExplicitHs(a.GetTotalNumHs() + a.GetNumRadicalElectrons())
                a.SetNumRadicalElectrons(0)
            mol.UpdatePropertyCache()
            clean_smiles = dm.to_smiles(mol)
            clean_mol = dm.to_mol(clean_smiles)
        except Exception:
            pass

        return clean_mol, clean_smiles

    def _is_valid(self, new_mol, new_smiles, need_substruct: Optional[dm.Mol] = None):
        """Check whether the mol or smiles is valid"""
        # not disconnected
        is_ok = new_smiles.count(".") == 0
        if not self.allow_cycle:
            is_ok &= not new_mol.HasSubstructMatch(self._ring_smarts)
        is_ok &= not (any(new_mol.HasSubstructMatch(x) for x in self._impossible_chemistry))
        if need_substruct is not None:
            is_ok &= new_mol.HasSubstructMatch(need_substruct)
        return is_ok

    def _has_correct_size(self, new_mol, expected_size):
        """Check whether the molecule has the correct size"""
        return new_mol.GetNumHeavyAtoms() == expected_size

    def __call__(self, *args, **kwargs):
        return self.enumerate(*args, **kwargs)

    def enumerate(
        self,
        mol,
        depth: Optional[int] = None,
        include_input: bool = True,
        protect_substruct: Optional[Chem.rdchem.Mol] = None,
        max_mols: Optional[int] = None,
    ):
        """Enumerate the list of isomers of the current mol

        Args:
            mol (Chem.rdchem.Mol): input molecule or smiles
            depth (int, optional): optional maximum depth of the breadth search (default=None)
            as_smiles (bool, optional): whether to return smiles strings of mols (default=False)
            include_input (bool, optional): whether to include the input molecule (default=True)
            protect_substruct (Chem.rdchem.Mol, optional): Optional substruct to protect
            max_mols (int, optional): maximum number of molecule to sample
        """

        mol = dm.to_mol(mol)
        if protect_substruct is not None:
            mol = dm.protect_atoms(mol, substruct=protect_substruct, in_place=False, as_smarts=True)
        smiles = dm.to_smiles(mol, isomeric=True, canonical=True)
        mol_size = mol.GetNumHeavyAtoms()
        source_mols = [(mol, 0)]
        seen = set([smiles])
        seen_inchi_key = set([dm.to_inchikey(smiles)])
        depth = depth or float("inf")
        max_mols = max_mols or float("inf")
        if include_input:
            yield smiles
        while len(source_mols) > 0 and max_mols > 0:
            curmol, curdepth = source_mols.pop(0)
            try:
                # we first try to kekulize the molecule
                Chem.Kekulize(curmol)
                curmol = Chem.AddHs(curmol, explicitOnly=False)
            except Exception as e:
                pass
            if curdepth >= depth:
                return
            with dm.without_rdkit_log():
                for k, rxn in self.rxn_cache.items():
                    if max_mols < 1:
                        break
                    for new_mols in rxn.RunReactants((curmol,)):
                        if max_mols < 1:
                            break
                        new_mol = new_mols[0]
                        # sanitize mol
                        tmp_sm = dm.to_smiles(new_mol)
                        new_mol, new_smiles = self._clean(new_mol)
                        # need to clean twice, not sure why
                        _, new_smiles = self._clean(new_mol)
                        new_smiles_inchikey = dm.to_inchikey(new_mol)

                        if (
                            new_mol is None
                            or new_smiles_inchikey in seen_inchi_key
                            or not self._has_correct_size(new_mol, mol_size)
                        ):
                            continue
                        if protect_substruct is not None:
                            new_mol = dm.protect_atoms(
                                new_mol,
                                substruct=protect_substruct,
                                in_place=False,
                            )
                        seen.add(new_smiles)
                        seen_inchi_key.add(new_smiles_inchikey)
                        source_mols.append((new_mol, curdepth + 1))
                        if self._is_valid(new_mol, new_smiles, protect_substruct):
                            yield new_smiles
                            max_mols -= 1
