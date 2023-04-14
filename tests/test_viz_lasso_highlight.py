import pytest
import datamol as dm


# The following tests are supposed to work and should not raise any errors
def test_original_working_solution_str():
    smi = "CO[C@@H](O)C1=C(O[C@H](F)Cl)C(C#N)=C1ONNC[NH3+]"
    smarts_list = "CONN"
    assert dm.lasso_highlight_image(smi, smarts_list)


# The following tests are supposed to work and should not raise any errors
def test_from_mol():
    smi = "CO[C@@H](O)C1=C(O[C@H](F)Cl)C(C#N)=C1ONNC[NH3+]"
    mol = dm.to_mol(smi)
    smarts_list = "CONN"
    assert dm.lasso_highlight_image(mol, smarts_list)


def test_original_working_solution_list_single_str():
    smi = "CO[C@@H](O)C1=C(O[C@H](F)Cl)C(C#N)=C1ONNC[NH3+]"
    smarts_list = ["CONN"]
    assert dm.lasso_highlight_image(smi, smarts_list)


def test_original_working_solution_list_str():
    smi = "CO[C@@H](O)C1=C(O[C@H](F)Cl)C(C#N)=C1ONNC[NH3+]"
    smarts_list = ["CONN", "N#CC~CO", "C=CON", "CONNCN"]
    assert dm.lasso_highlight_image(smi, smarts_list)


def test_original_working_solution_mol():
    smi = "CO[C@@H](O)C1=C(O[C@H](F)Cl)C(C#N)=C1ONNC[NH3+]"
    smarts_list = dm.to_mol("CONN")
    assert dm.lasso_highlight_image(smi, smarts_list)


def test_original_working_solution_list_single_mol():
    smi = "CO[C@@H](O)C1=C(O[C@H](F)Cl)C(C#N)=C1ONNC[NH3+]"
    smarts_list = [dm.to_mol("CONN")]
    assert dm.lasso_highlight_image(smi, smarts_list)


def test_original_working_solution_List_mol():
    smi = "CO[C@@H](O)C1=C(O[C@H](F)Cl)C(C#N)=C1ONNC[NH3+]"
    smarts_list = [dm.to_mol("CONN"), dm.to_mol("N#CC~CO"), dm.to_mol("C=CON"), dm.to_mol("CONNCN")]
    assert dm.lasso_highlight_image(smi, smarts_list)


def test_wokring_solution_with_more_structures_than_colors():
    smi = "CO[C@@H](O)C1=C(O[C@H](F)Cl)C(C#N)=C1ONNC[NH3+]"
    smarts_list = ["CONN", "N#CC~CO", "C=CON", "CONNCN", "FCCl", "OCO", "N#C", "N#CC", "CC#N"]
    assert dm.lasso_highlight_image(smi, smarts_list)


def test_drawing_options():
    smi = "CO[C@@H](O)C1=C(O[C@H](F)Cl)C(C#N)=C1ONNC[NH3+]"
    smarts_list = ["CONN", "N#CC~CO", "C=CON", "CONNCN", "FCCl", "OCO", "N#C", "N#CC", "CC#N"]
    assert dm.lasso_highlight_image(smi, smarts_list, bondLineWidth=15)


def test_wrong_drawing_options():
    smi = "CO[C@@H](O)C1=C(O[C@H](F)Cl)C(C#N)=C1ONNC[NH3+]"
    smarts_list = ["CONN", "N#CC~CO", "C=CON", "CONNCN", "FCCl", "OCO", "N#C", "N#CC", "CC#N"]

    with pytest.raises(ValueError):
        dm.lasso_highlight_image(smi, smarts_list, bondLineWidthXXXXXXX=15)


def test_input_mol_is_none():
    smarts_list = ["CONN", "N#CC~CO", "C=CON", "CONNCN", "FCCl", "OCO", "N#C", "N#CC", "CC#N"]

    with pytest.raises(ValueError):
        dm.lasso_highlight_image(None, smarts_list)


# The following tests are supposed to raise errors and will be check to ensure they do
def test_canvas_input_error():
    with pytest.raises(ValueError):
        smi = "CO[C@@H](O)C1=C(O[C@H](F)Cl)C(C#N)=C1ONNC[NH3+]"
        smarts_list = ["CONN", "N#CC~CO", "C=CON", "CONNCN"]
        dm.lasso_highlight_image(smi, smarts_list, (5001, 400))


def test_search_input_error_empty_str():
    with pytest.raises(ValueError):
        smi = "CO[C@@H](O)C1=C(O[C@H](F)Cl)C(C#N)=C1ONNC[NH3+]"
        smarts_list = ""
        dm.lasso_highlight_image(smi, smarts_list)


def test_search_input_error_empty_list():
    # should still go through but just print out the structure without any highlights
    smi = "CO[C@@H](O)C1=C(O[C@H](F)Cl)C(C#N)=C1ONNC[NH3+]"
    smarts_list = []
    assert dm.lasso_highlight_image(smi, smarts_list)


def test_search_input_error_None():
    with pytest.raises(ValueError):
        smi = "CO[C@@H](O)C1=C(O[C@H](F)Cl)C(C#N)=C1ONNC[NH3+]"
        smarts_list = None
        dm.lasso_highlight_image(smi, smarts_list)


def test_target_input_error_empty_str():
    with pytest.raises(ValueError):
        smi = ""
        smarts_list = ["CONN", "N#CC~CO", "C=CON", "CONNCN"]
        dm.lasso_highlight_image(smi, smarts_list)


def test_target_input_error_None():
    with pytest.raises(ValueError):
        smi = None
        smarts_list = ["CONN", "N#CC~CO", "C=CON", "CONNCN"]
        dm.lasso_highlight_image(smi, smarts_list)


def test_search_input_error_smarts_no_substructure():
    # This test should still continue but will just print out a structure without any highlights and a warning
    smi = "CO[C@@H](O)C1=C(O[C@H](F)Cl)C(C#N)=C1ONNC[NH3+]"
    smarts_list = ["CCCCCC"]
    assert dm.lasso_highlight_image(smi, smarts_list)


# testing using <class 'IPython.core.display.SVG'>" == str(type(img)) so to not bring in IPython
# as a dependency for the tests
def test_SVG_is_returned_explicit():
    smi = "CO[C@@H](O)C1=C(O[C@H](F)Cl)C(C#N)=C1ONNC[NH3+]"
    smarts_list = ["CC"]
    img = dm.lasso_highlight_image(smi, smarts_list, use_svg=True)
    assert isinstance(img, str)


def test_SVG_is_returned_implict():
    smi = "CO[C@@H](O)C1=C(O[C@H](F)Cl)C(C#N)=C1ONNC[NH3+]"
    smarts_list = ["CC"]
    img = dm.lasso_highlight_image(smi, smarts_list)
    assert isinstance(img, str)


def test_PNG_is_returned():
    smi = "CO[C@@H](O)C1=C(O[C@H](F)Cl)C(C#N)=C1ONNC[NH3+]"
    smarts_list = ["CC"]
    img = dm.lasso_highlight_image(smi, smarts_list, use_svg=False)

    from PIL import Image

    assert isinstance(img, Image.Image)
