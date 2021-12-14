import exago
from test_utilities import *
from check_preconditions import check_preconditions
check_preconditions()


@exago_test
def test_initialize():
    '''Test calling ExaGOInitialize'''
    exago.initialize("opflow")


@exago_test
def test_creating_opflow():
    '''Testing creation of opflow object'''
    opf = exago.opf()
    print(opf)


@exago_test
def test_tolerance_opflow():
    '''Testing setting tolerance in opflow object'''
    opf = exago.opf()
    opf.set_tolerance(1e-6)
    tol = opf.get_tolerance()
    assert tol == 1e-6
    opf.set_tolerance(1e-5)
    tol = opf.get_tolerance()
    assert tol == 1e-5


@exago_test
def test_set_loadloss_penalty_opflow():
    '''Testing setting load loss penalty in opflow object'''
    opf = exago.opf()
    opf.set_loadloss_penalty(999)


@exago_test
def test_set_powerimbalance_penalty_opflow():
    '''Testing setting power imbalance penalty in opflow object'''
    opf = exago.opf()
    opf.set_powerimbalance_penalty(999)


@exago_test
def test_solve_opflow():
    '''Testing opflow solve'''
    opf = exago.opf()
    print(opf)
    path = exago.get_datafile_path()
    opf.read_mat_power_data(os.path.join(path, 'case9', 'case9mod.m'))
    opf.solve()


@exago_test
def test_get_datafile_path():
    '''Test retrieving datafile path'''
    path = exago.get_datafile_path()
    assert isinstance(path, str)


@exago_test
def test_finalize():
    '''Test calling ExaGOFinalize'''
    exago.finalize()


exago_run_all_tests(__file__)
