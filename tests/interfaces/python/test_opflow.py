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
    opf.setup_ps()
    opf.ps_set_gen_power_limits(3, "1 ", 65, 0, exago_ignore, exago_ignore)
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
def test_set_model():
    '''Test setting opflow model'''
    opf = OPFLOW()
    opf.dont_finalize()
    opf.set_model('POWER_BALANCE_POLAR')

@exago_test
def test_set_solver():
    '''Test setting opflow solver'''
    opf = OPFLOW()
    opf.dont_finalize()
    opf.read_mat_power_data(
        os.path.join(config.prefix(), 'share', 'exago', 'datafiles', 'case9', 'case9mod.m'))
    opf.set_model('POWER_BALANCE_CARTESIAN')
    opf.set_solver('TAO')
    opf.solve()

@exago_test
def test_set_initialization_type():
    '''Test setting opflow intialization type'''
    opf = OPFLOW()
    opf.dont_finalize()
    opf.read_mat_power_data(
        os.path.join(config.prefix(), 'share', 'exago', 'datafiles', 'case9', 'case9mod.m'))
    opf.set_initialization('MIDPOINT')
    opf.solve()

@exago_test
def test_set_ignore_lineflow_constraints():
    '''Test setting opflow ignore_lineflow_constraints'''
    opf = OPFLOW()
    opf.dont_finalize()
    opf.read_mat_power_data(
        os.path.join(config.prefix(), 'share', 'exago', 'datafiles', 'case9', 'case9mod.m'))
    opf.set_ignore_lineflow_constraints(True)
    opf.solve()

@exago_test
def test_set_include_loadloss():
    '''Test setting opflow ignore_lineflow_constraints'''
    opf = OPFLOW()
    opf.dont_finalize()
    opf.read_mat_power_data(
        os.path.join(config.prefix(), 'share', 'exago', 'datafiles', 'case9', 'case9mod.m'))
    opf.set_include_loadloss(True)
    opf.solve()

@exago_test
def test_set_include_powerimbalance():
    '''Test setting opflow include powerimbalance variables'''
    opf = OPFLOW()
    opf.dont_finalize()
    opf.read_mat_power_data(
        os.path.join(config.prefix(), 'share', 'exago', 'datafiles', 'case9', 'case9mod.m'))
    opf.set_include_powerimbalance(True)
    opf.solve()

@exago_test
def test_set_loadloss_penalty():
    '''Test setting opflow loadloss penalty'''
    opf = OPFLOW()
    opf.dont_finalize()
    opf.read_mat_power_data(
        os.path.join(config.prefix(), 'share', 'exago', 'datafiles', 'case9', 'case9mod.m'))
    opf.set_include_loadloss(True)
    opf.set_loadloss_penalty(999)
    opf.set_initialization('FROMFILE')
    opf.set_genbusvoltage('VARIABLE_WITHIN_BOUNDS')
    opf.solve()

@exago_test
def test_set_include_powerimbalance():
    '''Test setting opflow include powerimbalance variables'''
    opf = OPFLOW()
    opf.dont_finalize()
    opf.read_mat_power_data(
        os.path.join(config.prefix(), 'share', 'exago', 'datafiles', 'case9', 'case9mod.m'))
    opf.set_include_powerimbalance(True)
    opf.set_include_loadloss(True)
    opf.set_loadloss_penalty(999)
    opf.set_powerimbalance_penalty(999)
    opf.set_initialization('FROMFILE')
    opf.set_genbusvoltage('VARIABLE_WITHIN_BOUNDS')
    opf.solve()

@exago_test
def test_set_genbusvoltage():
    '''Test setting opflow genbusvoltage mode'''
    opf = OPFLOW()
    opf.dont_finalize()
    opf.read_mat_power_data(
        os.path.join(config.prefix(), 'share', 'exago', 'datafiles', 'case9', 'case9mod.m'))
    opf.set_genbusvoltage('VARIABLE_WITHIN_BOUNDS')
    opf.solve()

@exago_test
def test_set_has_gensetpoint():
    '''Test setting opflow has gensetpoint variable'''
    opf = OPFLOW()
    opf.dont_finalize()
    opf.read_mat_power_data(
        os.path.join(config.prefix(), 'share', 'exago', 'datafiles', 'case9', 'case9mod.m'))
    opf.set_has_gensetpoint(True)
    opf.solve()

@exago_test
def test_set_objective_function():
    '''Test setting opflow objective function'''
    opf = OPFLOW()
    opf.dont_finalize()
    opf.read_mat_power_data(
        os.path.join(config.prefix(), 'share', 'exago', 'datafiles', 'case9', 'case9mod.m'))
    opf.set_objective_function('MIN_GENSETPOINT_DEVIATION')
    opf.solve()

@exago_test
def test_set_use_agc():
    '''Test setting opflow AGC for generator real power dispatch'''
    opf = OPFLOW()
    opf.dont_finalize()
    opf.read_mat_power_data(
        os.path.join(config.prefix(), 'share', 'exago', 'datafiles', 'case9', 'case9mod.m'))
    opf.set_use_agc(True)
    opf.solve()

@exago_test
def test_set_tolerance():
    '''Test setting tolerance for optimization solver'''
    opf = OPFLOW()
    opf.dont_finalize()
    opf.read_mat_power_data(
        os.path.join(config.prefix(), 'share', 'exago', 'datafiles', 'case9', 'case9mod.m'))
    opf.set_tolerance(0.000001)
    opf.solve()

@exago_test
def test_set_hiop_compute_mode():
    '''Test setting hiop compute mode'''
    opf = OPFLOW()
    opf.dont_finalize()
    opf.read_mat_power_data(
        os.path.join(config.prefix(), 'share', 'exago', 'datafiles', 'case9', 'case9mod.m'))
    opf.set_model('POWER_BALANCE_HIOP')
    opf.set_solver('HIOP')
    opf.set_hiop_compute_mode('cpu')
    opf.set_genbusvoltage('VARIABLE_WITHIN_BOUNDS')
    opf.solve()

@exago_test
def test_set_hiop_verbosity():
    '''Test setting hiop verbosity level'''
    opf = OPFLOW()
    opf.dont_finalize()
    opf.read_mat_power_data(
        os.path.join(config.prefix(), 'share', 'exago', 'datafiles', 'case9', 'case9mod.m'))
    opf.set_model('POWER_BALANCE_HIOP')
    opf.set_solver('HIOP')
    opf.set_hiop_compute_mode('cpu')
    opf.set_genbusvoltage('VARIABLE_WITHIN_BOUNDS')
    opf.set_hiop_verbosity(3)
    opf.solve()

@exago_test
def test_print_save_output():
    '''Testing printing and saving output'''
    opf = OPFLOW()
    opf.dont_finalize()
    opf.read_mat_power_data(
        os.path.join(config.prefix(), 'share', 'exago', 'datafiles', 'case9', 'case9mod.m'))
    opf.solve()
    opf.print_solution()
    filepath = os.path.join(config.prefix(), '..', 'tmp.output')
    opf.save_solution('CSV', filepath) 

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
