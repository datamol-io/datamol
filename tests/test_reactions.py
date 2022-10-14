import pytest
import datamol as dm


def test_convert_attach_to_isotope():
    smiles = "O=C(NC[*:1])Cc1ccnc(OCc2ccccc2)c1"

    isotope = dm.reactions.convert_attach_to_isotope(mol_or_smiles=smiles)

    assert dm.to_smiles(isotope, canonical=True) == dm.to_smiles(
        dm.to_mol("O=C(Cc1ccnc(OCc2ccccc2)c1)NC[1*]"), canonical=True
    )


def test_num_attachment_points():
    smiles = "O=C(NC[*:1])Cc1c([*:2])nc(OCc2ccccc2)c1"
    assert dm.reactions.num_attachment_points(smiles) == 2


def test_open_attach_points():
    mol = dm.to_mol("CC")
    smiles_open = dm.to_smiles(
        dm.reactions._attachments.open_attach_points(mol, bond_type=dm.SINGLE_BOND), canonical=True
    )
    assert smiles_open == dm.to_smiles(dm.to_mol("*CC*"), canonical=True)


def test_rxn_from_smarts():
    smarts = "[1*][*:1].[1*][*:2]>>[*:1][*:2]"
    rxn = dm.reactions.rxn_from_smarts(smarts)
    assert len(list(rxn.GetReactants())) == 2


def test_rxn_from_block():
    block = """$RXN

      ISIS     090220091539

  2  1
$MOL

  -ISIS-  09020915392D

  2  1  1  0  0  0  0  0  0  0999 V2000
   -2.0744    0.1939    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5440   -0.1592    0.0000 R#  0  0  0  0  0  0  0  0  0  1  0  0
  1  2  1  0  0  0  0
  1 F    2  17  35
V    1 halogen
M  RGP  1   2   1
M  ALS   1  2 F Cl  Br
M  END
$MOL

  -ISIS-  09020915392D

  2  1  0  0  0  0  0  0  0  0999 V2000
    2.8375   -0.2500    0.0000 R#  0  0  0  0  0  0  0  0  0  2  0  0
    3.3463    0.0438    0.0000 N   0  0  0  0  0  0  0  0  0  3  0  0
  1  2  1  0  0  0  0
V    2 amine.primary
M  RGP  1   1   2
M  END
$MOL

  -ISIS-  09020915392D

  3  2  0  0  0  0  0  0  0  0999 V2000
   13.5792    0.0292    0.0000 N   0  0  0  0  0  0  0  0  0  3  0  0
   14.0880    0.3229    0.0000 R#  0  0  0  0  0  0  0  0  0  1  0  0
   13.0704    0.3229    0.0000 R#  0  0  0  0  0  0  0  0  0  2  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
M  RGP  2   2   1   3   2
M  END"""
    rxn = dm.reactions.rxn_from_block(block)
    assert len(list(rxn.GetReactants())) == 2


def test_is_reaction_ok():
    smarts = "[1*][*:1].[1*][*:2]>>[*:1][*:2]"
    rxn = dm.reactions.rxn_from_smarts(smarts)
    assert dm.reactions.is_reaction_ok(rxn) == True


def test_sanitize():
    block = """$RXN

      ISIS     090220091539

  2  1
$MOL

  -ISIS-  09020915392D

  2  1  1  0  0  0  0  0  0  0999 V2000
   -2.0744    0.1939    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5440   -0.1592    0.0000 R#  0  0  0  0  0  0  0  0  0  1  0  0
  1  2  1  0  0  0  0
  1 F    2  17  35
V    1 halogen
M  RGP  1   2   1
M  ALS   1  2 F Cl  Br
M  END
$MOL

  -ISIS-  09020915392D

  2  1  0  0  0  0  0  0  0  0999 V2000
    2.8375   -0.2500    0.0000 R#  0  0  0  0  0  0  0  0  0  2  0  0
    3.3463    0.0438    0.0000 Na   0  0  0  0  0  0  0  0  0  3  0  0
  1  2  1  0  0  0  0
V    2 amine.primary
M  RGP  1   1   2
M  END
$MOL

  -ISIS-  09020915392D

  3  2  0  0  0  0  0  0  0  0999 V2000
   13.5792    0.0292    0.0000 Na   0  0  0  0  0  0  0  0  0  3  0  0
   14.0880    0.3229    0.0000 R#  0  0  0  0  0  0  0  0  0  1  0  0
   13.0704    0.3229    0.0000 R#  0  0  0  0  0  0  0  0  0  2  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
M  RGP  2   2   1   3   2
M  END"""
    with pytest.raises(ValueError):
        dm.reactions.rxn_from_block(rxn_block_or_file=block, sanitize=True)


def test_apply_reaction():
    # 1 attach point
    rxn = dm.reactions.rxn_from_smarts("[1*][*:1].[1*][*:2]>>[*:1][*:2]")
    frag = dm.to_mol("CC(O[1*])C(F)(F)F")
    scf = dm.to_mol("O=C(NC[1*])NCc1ccnc(OCc2ccccc2)c1")
    prod = dm.reactions.apply_reaction(rxn=rxn, reactants=(scf, frag), single_output=True)
    assert len([prod]) == 1
    assert isinstance(prod, dm.Mol) == True

    prod = dm.reactions.apply_reaction(
        rxn=rxn, reactants=(scf, frag), single_output=False, as_smiles=True
    )
    assert len(prod) >= 1
    assert isinstance(prod, list) == True


def test_can_react():
    rxn = dm.reactions.rxn_from_smarts(
        "[*:1]O[C:2](=[O])[#6:3].[Nh:4]>>[N:4][C:2](=[O])[#6:3].[*:1][Oh]"
    )
    mol = dm.to_mol("Nc1ccccc1")
    assert dm.reactions.can_react(rxn, mol) == True

    mol = dm.to_mol("CC")
    assert dm.reactions.can_react(rxn, mol) == False


def test_inverse_reaction():
    rxn = dm.reactions.rxn_from_smarts("[1*][*:1].[1*][*:2]>>[*:1][*:2]")
    rxn_r = dm.reactions.inverse_reaction(rxn)
    assert len(list(rxn_r.GetReactants())) == 1
    assert len(list(rxn_r.GetProducts())) == 2
