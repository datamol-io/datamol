import pytest

import datamol as dm


@dm.no_rdkit_log
def no_log_to_mol(smiles):
    return dm.to_mol(smiles)


def check_logs_are_shown(capfd):
    smiles = "fake_smiles"
    dm.to_mol(smiles)
    _, err = capfd.readouterr()
    assert "SMILES Parse Error" in err


def check_logs_are_not_shown(capfd):
    smiles = "fake_smiles"
    dm.to_mol(smiles)
    _, err = capfd.readouterr()
    assert err == ""


def check_logs_are_not_shown_deco(capfd):
    smiles = "fake_smiles"
    no_log_to_mol(smiles)
    _, err = capfd.readouterr()
    assert err == ""


@pytest.mark.skip_platform("win")
def test_rdkit_log(capfd):
    """Test multiple rdkit log scenarios."""

    check_logs_are_shown(capfd)
    check_logs_are_not_shown_deco(capfd)

    check_logs_are_shown(capfd)
    with dm.without_rdkit_log():
        check_logs_are_not_shown(capfd)
    check_logs_are_shown(capfd)

    dm.disable_rdkit_log()
    check_logs_are_not_shown(capfd)

    dm.enable_rdkit_log()
    check_logs_are_shown(capfd)

    dm.disable_rdkit_log()
    with dm.without_rdkit_log():
        check_logs_are_not_shown(capfd)
    check_logs_are_not_shown(capfd)


@pytest.mark.skip_platform("win")
def test_rdkit_log_enable(capfd):
    dm.enable_rdkit_log()

    with dm.without_rdkit_log():
        check_logs_are_not_shown(capfd)

    with dm.without_rdkit_log(enable=False):
        check_logs_are_shown(capfd)

    check_logs_are_shown(capfd)
