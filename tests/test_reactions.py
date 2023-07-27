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
    assert dm.reactions.is_reaction_ok(rxn) is True


def test_is_reaction_ok_with_logs(caplog):
    smarts = "[1*][*:1].[1*][*:2]>>[*:1][*:2]"
    rxn = dm.reactions.rxn_from_smarts(smarts)
    assert dm.reactions.is_reaction_ok(rxn, enable_logs=True) is True

    assert "Number of" in caplog.text


def test_apply_reaction():
    # 1 attach point
    rxn = dm.reactions.rxn_from_smarts("[1*][*:1].[1*][*:2]>>[*:1][*:2]")
    frag = dm.to_mol("CC(O[1*])C(F)(F)F")
    scf = dm.to_mol("O=C(NC[1*])NCc1ccnc(OCc2ccccc2)c1")
    prod = dm.reactions.apply_reaction(
        rxn=rxn, reactants=(scf, frag), product_index=0, single_product_group=True
    )
    assert len([prod]) == 1
    assert isinstance(prod, dm.Mol) is True

    prod = dm.reactions.apply_reaction(
        rxn=rxn, reactants=(scf, frag), single_product_group=False, as_smiles=True
    )
    assert len(prod) >= 1
    assert isinstance(prod, list) is True


def test_apply_reaction_invalid():
    # 1 attach point
    rxn = dm.reactions.rxn_from_smarts("[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]")
    frag = dm.to_mol("CC(O)C(F)(F)F")
    scf = dm.to_mol("O=C(NC)NCc1ccnc(OCc2cc[2*]ccc2)c1")
    prod = dm.reactions.apply_reaction(
        rxn=rxn, reactants=(scf, frag), product_index=0, single_product_group=True
    )
    assert isinstance(prod, list)
    assert len(prod) == 0


def test_can_react():
    rxn = dm.reactions.rxn_from_smarts(
        "[*:1]O[C:2](=[O])[#6:3].[Nh:4]>>[N:4][C:2](=[O])[#6:3].[*:1][Oh]"
    )
    mol = dm.to_mol("Nc1ccccc1")
    assert dm.reactions.can_react(rxn, mol) is True

    mol = dm.to_mol("CC")
    assert dm.reactions.can_react(rxn, mol) is False


def test_inverse_reaction():
    rxn = dm.reactions.rxn_from_smarts("[1*][*:1].[1*][*:2]>>[*:1][*:2]")
    rxn_r = dm.reactions.inverse_reaction(rxn)
    assert len(list(rxn_r.GetReactants())) == 1
    assert len(list(rxn_r.GetProducts())) == 2


def test_select_reaction_output():
    smiles = (
        (
            "Cc1cnc(CNc(cccc2-c3cn(CC(C4)CC4O)c4ncnc(N)c34)c2F)s1",
            "Cc1cnc(CNc(cccc2-c3cn(C(C4)CC4O)c4ncnc(N)c34)c2F)s1",
        ),
        (
            "Nc1c(c(-c2cccc(NCc3nccs3)c2F)cn2Cc3cn(C4CC5(COCC5)OCC4)nn3)c2ncn1",
            "C[C@@H](c1cn(CC2CCCCC2)nn1)N(Cc1nccs1)c(cccc1-c2c[nH]c3ncnc(N)c23)c1F",
        ),
        (
            "Nc1c(c(-c2cccc(NCc3nccs3)c2F)cn2Cc3cn(CCN4CCOCC4)nn3)c2ncn1",
            "Cc1cnc(C(Nc(cccc2-c3cn(C[C@H](C4)C[C@@H]4O)c4ncnc(N)c34)c2F)=O)s1",
        ),
    )
    product = tuple(tuple([dm.to_mol(p[0]), dm.to_mol(p[1])]) for p in smiles)
    prod = dm.reactions.select_reaction_output(
        product,
        product_index=None,
        single_product_group=False,
        rm_attach=True,
        as_smiles=False,
        sanitize=True,
    )
    assert len(prod) == len(smiles)

    product = tuple(tuple([dm.to_mol(p[0]), dm.to_mol(p[1])]) for p in smiles)
    prod = dm.reactions.select_reaction_output(
        product,
        product_index=[0, 1],
        single_product_group=False,
        rm_attach=True,
        as_smiles=False,
        sanitize=True,
    )
    assert len(prod) == len(smiles)

    product = tuple(tuple([dm.to_mol(p[0]), dm.to_mol(p[1])]) for p in smiles)
    prod = dm.reactions.select_reaction_output(
        product,
        product_index=1,
        single_product_group=True,
        rm_attach=True,
        as_smiles=False,
        sanitize=True,
    )
    assert isinstance(prod, dm.mol.Mol) is True


def test_select_reaction_output_invalid():
    smiles = (
        (
            "Cc1cnc(CNc(cccc2-c3cn(CC(C4)CC4O)c4ncnc(N)c34)c2F)s1",
            "ABCDEFG",  # invalid mol
        ),
        (
            "HIJKLMN",  # invalid mol
            "C[C@@H](c1cn(CC2CCCCC2)nn1)N(Cc1nccs1)c(cccc1-c2c[nH]c3ncnc(N)c23)c1F",
        ),
        (
            "Nc1c(c(-c2cccc(NCc3nccs3)c2F)cn2Cc3cn(CCN4CCOCC4)nn3)c2ncn1",
            "Cc1cnc(C(Nc(cccc2-c3cn(C[C@H](C4)C[C@@H]4O)c4ncnc(N)c34)c2F)=O)s1",
        ),
    )

    product = tuple(tuple([dm.to_mol(p[0]), dm.to_mol(p[1])]) for p in smiles)
    prod = dm.reactions.select_reaction_output(
        product=product,
        product_index=1,
        single_product_group=False,
        rm_attach=True,
        as_smiles=True,
        sanitize=True,
    )
    assert prod[0] is None
    assert isinstance(prod[1], str)
    assert isinstance(prod[2], str)
