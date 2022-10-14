import datamol as dm

REACTION_BLOCK = """$RXN

      RDKit

  2  1
$MOL

     RDKit          2D

  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 *   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    0.7500    0.0000 *   0  0  0  0  0  0  0  0  0  1  0  0
  1  2  6  0
V    1 [1*]
V    2 [*:1]
M  END
$MOL

     RDKit          2D

  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 *   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    0.7500    0.0000 *   0  0  0  0  0  0  0  0  0  2  0  0
  1  2  6  0
V    1 [1*]
V    2 [*:2]
M  END
$MOL

     RDKit          2D

  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 *   0  0  0  0  0  0  0  0  0  1  0  0
    1.2990    0.7500    0.0000 *   0  0  0  0  0  0  0  0  0  2  0  0
  1  2  6  0
V    1 [*:1]
V    2 [*:2]
M  END"""


def test_convert_attach_to_isotope():
    smiles = "O=C(NC[*:1])Cc1ccnc(OCc2ccccc2)c1"

    isotope = dm.reactions.convert_attach_to_isotope(mol_or_smiles=smiles)

    assert dm.to_smiles(isotope, canonical=True) == dm.to_smiles(
        dm.to_mol("O=C(Cc1ccnc(OCc2ccccc2)c1)NC[1*]"), canonical=True
    )

    # same isotope
    dm.reactions.convert_attach_to_isotope(mol_or_smiles=smiles, same_isotope=True)

    isotope2 = dm.reactions.convert_attach_to_isotope(mol_or_smiles=smiles, as_smiles=True)
    assert isotope2 == "O=C(Cc1ccnc(OCc2ccccc2)c1)NC[1*]"


def test_num_attachment_points():
    smiles = "O=C(NC[*:1])Cc1c([*:2])nc(OCc2ccccc2)c1"
    assert dm.reactions.num_attachment_points(smiles) == 2

    assert dm.reactions.num_attachment_points(dm.from_smarts(smiles)) == 2


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


def test_rxn_to_smarts():
    smarts = "[1*][*:1].[1*][*:2]>>[*:1][*:2]"
    rxn = dm.reactions.rxn_from_smarts(smarts)

    smarts2 = dm.reactions.rxn_to_smarts(rxn)
    assert smarts == smarts2


def test_rxn_from_block():
    rxn = dm.reactions.rxn_from_block(REACTION_BLOCK)
    assert rxn.GetNumProductTemplates() == 1
    assert rxn.GetNumReactantTemplates() == 2


def test_rxn_from_block_file(tmp_path):
    block_path = tmp_path / "block.block"

    with open(block_path, "w") as f:
        f.write(REACTION_BLOCK)

    rxn = dm.reactions.rxn_from_block_file(block_path)
    assert rxn.GetNumProductTemplates() == 1
    assert rxn.GetNumReactantTemplates() == 2


def test_rxn_to_block():
    rxn = dm.reactions.rxn_from_smarts("[*:1][*:2]>>[*:1].[*:2]")
    block = dm.reactions.rxn_to_block(rxn)
    rxn2 = dm.reactions.rxn_from_block(block)
    assert dm.reactions.rxn_to_smarts(rxn2) == dm.reactions.rxn_to_smarts(rxn)


def test_rxn_to_block_file(tmp_path):
    block_path = tmp_path / "block.block"
    rxn = dm.reactions.rxn_from_smarts("[*:1][*:2]>>[*:1].[*:2]")
    dm.reactions.rxn_to_block_file(rxn, block_path)

    rxn2 = dm.reactions.rxn_from_block_file(block_path)
    assert dm.reactions.rxn_to_smarts(rxn2) == dm.reactions.rxn_to_smarts(rxn)


def test_is_reaction_ok():
    smarts = "[1*][*:1].[1*][*:2]>>[*:1][*:2]"
    rxn = dm.reactions.rxn_from_smarts(smarts)
    assert dm.reactions.is_reaction_ok(rxn) == True


def test_is_reaction_ok_with_logs(caplog):
    smarts = "[1*][*:1].[1*][*:2]>>[*:1][*:2]"
    rxn = dm.reactions.rxn_from_smarts(smarts)
    assert dm.reactions.is_reaction_ok(rxn, enable_logs=True) == True

    assert "Number of" in caplog.text


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
