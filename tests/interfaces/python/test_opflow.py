from exago.opflow import OPFLOW
from exago import config
from test_utilities import *
from check_preconditions import check_preconditions
from ctypes import *
check_preconditions()

# This is a workaround for not having access to EXAGO_IGNORE
exago_ignore = -1000000

@exago_test
def test_opflow_create():
    '''Test creation of powerflow object'''
    opf = OPFLOW()
    opf.dont_finalize()

@exago_test
def test_opflow_read_mat_power_data():
    '''Test powerflow reading data in mat power file format'''
    opf = OPFLOW()
    opf.dont_finalize()
    opf.read_mat_power_data(
        os.path.join(config.prefix(), 'share', 'exago', 'datafiles', 'case9', 'case9mod.m'))

@exago_test
def test_opflow_solve():
    '''Test solving a powerflow problem'''
    opf = OPFLOW()
    opf.dont_finalize()
    opf.read_mat_power_data(
        os.path.join(config.prefix(), 'share', 'exago', 'datafiles', 'case9', 'case9mod.m'))
    opf.solve()

@exago_test
def test_solution_to_ps():
    '''Test storing opflow solution to PS structure internally'''
    opf = OPFLOW()
    opf.dont_finalize()
    opf.read_mat_power_data(
        os.path.join(config.prefix(), 'share', 'exago', 'datafiles', 'case9', 'case9mod.m'))
    opf.solve()
    opf.solution_to_ps()

@exago_test
def test_objective_function():
    '''Test getting the objective function from the opflow solution'''
    opf = OPFLOW()
    opf.dont_finalize()
    opf.read_mat_power_data(
        os.path.join(config.prefix(), 'share', 'exago', 'datafiles', 'case9', 'case9mod.m'))
    opf.solve()
    opf.solution_to_ps()
    objective = opf.objective_function

@exago_test
def test_setup_ps():
    '''Test setting up PS structure internal to opflow object'''
    opf = OPFLOW()
    opf.dont_finalize()
    opf.read_mat_power_data(
        os.path.join(config.prefix(), 'share', 'exago', 'datafiles', 'case9', 'case9mod.m'))
    opf.setup_ps()

@exago_test
def test_ps_set_gen_power_limits():
    '''Test setting PS structure through wrapper function'''
    opf = OPFLOW()
    opf.dont_finalize()
    opf.read_mat_power_data(
        os.path.join(config.prefix(), 'share', 'exago', 'datafiles', 'case9', 'case9mod.m'))
    opf.setup_ps()
    opf.ps_set_gen_power_limits(3, "1 ", 65, 0, exago_ignore, exago_ignore)
    opf.solve()

@exago_test
def test_get_gen_dispatch():
    '''Test getting the generator dispatch after opflow solve'''
    opf = OPFLOW()
    opf.dont_finalize()
    opf.read_mat_power_data(
        os.path.join(config.prefix(), 'share', 'exago', 'datafiles', 'case9', 'case9mod.m'))
    opf.setup_ps()
    opf.solve()
    opf.solution_to_ps()
    pg1set, qg1set = opf.get_gen_dispatch(1, "1 ")
    pg2set, qg2set = opf.get_gen_dispatch(2, "1 ")

@exago_test
def test_opf2_driver():
    '''Test opflow2 driver'''
    opf = OPFLOW()
    # We will finalize for second opf automatically
    # opf.dont_finalize()
    opf.read_mat_power_data(
        os.path.join(config.prefix(), 'share', 'exago', 'datafiles', 'case9', 'case9mod.m'))
    opf.setup_ps()
    opf.ps_set_gen_power_limits(3, "1 ", 65, 0, exago_ignore, exago_ignore)
    opf.solve()
    opf.solution_to_ps()
    obj1 = opf.objective_function
    pg1set, qg1set = opf.get_gen_dispatch(1, "1 ")
    pg2set, qg2set = opf.get_gen_dispatch(2, "1 ")
    opf2 = OPFLOW()
    opf2.read_mat_power_data(
        os.path.join(config.prefix(), 'share', 'exago', 'datafiles', 'case9', 'case9mod.m'))
    opf2.setup_ps()
    opf2.ps_set_gen_power_limits(3, "1 ", 105, 0, exago_ignore, exago_ignore)
    opf2.ps_set_gen_power_limits(1, "1 ", pg1set, pg1set, exago_ignore, exago_ignore)
    opf2.ps_set_gen_power_limits(2, "1 ", pg2set, pg2set, exago_ignore, exago_ignore)
    opf2.solve()
    opf2.solution_to_ps()
    obj2 = opf2.objective_function
    
exago_run_all_tests(__file__)
