import exago
from test_utilities import *
from check_preconditions import check_preconditions
check_preconditions()

exago_ignore = -1000000


@exago_test
def test_initialize():
    '''Test calling ExaGOInitialize'''
    exago.initialize("opflow")


@exago_test
def test_creating_opflow():
    '''Testing creation of opflow object'''
    opf = exago.OPFLOW()


@exago_test
def test_get_prefix():
    '''Test retrieving datafile path'''
    path = exago.prefix()
    assert isinstance(path, str)


@exago_test
def test_set_obj_func_ints():
    '''Testing setting objective function via integers'''
    opf = exago.OPFLOW()
    types = opf.get_objective_types()
    # List of possible enum values (object list)
    print(types)
    for i in range(0, len(types)):
        opf.set_objective_type(types[i])
        obj = opf.get_objective_type()
        assert isinstance(obj, exago.OPFLOWObjectiveType)
        assert obj == types[i]


@exago_test
def test_set_obj_func_strings():
    '''Testing setting objective function via strings'''
    opf = exago.OPFLOW()
    types = list(exago.OPFLOWObjectiveType.__members__.items())
    # List of tuple pairs of stringy values and corresponding enum value
    print(types)
    for i in range(0, len(types)):
        opf.set_objective_type(types[i][0])
        obj = opf.get_objective_type()
        assert isinstance(obj, exago.OPFLOWObjectiveType)
        assert obj == types[i][1]


@exago_test
def test_set_initialization_ints():
    '''Testing setting initialization type via integers'''
    opf = exago.OPFLOW()
    types = opf.get_initialization_types()
    # List of possible enum values (object list)
    print(types)
    for i in range(0, len(types)):
        opf.set_initialization_type(types[i])
        ini = opf.get_initialization_type()
        assert isinstance(ini, exago.OPFLOWInitializationType)
        assert ini == types[i]


@exago_test
def test_set_initialization_strings():
    '''Testing setting initialization type via strings'''
    opf = exago.OPFLOW()
    types = list(exago.OPFLOWInitializationType.__members__.items())
    # List of tuple pairs of stringy values and corresponding enum value
    print(types)
    for i in range(0, len(types)):
        opf.set_initialization_type(types[i][0])
        ini = opf.get_initialization_type()
        assert isinstance(ini, exago.OPFLOWInitializationType)
        assert ini == types[i][1]


@exago_test
def test_set_genbusvoltage_ints():
    '''Testing setting gen bus voltage type via ints'''
    opf = exago.OPFLOW()
    types = opf.get_gen_bus_voltage_types()
    # List of possible enum values (object list)
    print(types)
    for i in range(0, len(types)):
        opf.set_gen_bus_voltage_type(types[i])
        volt = opf.get_gen_bus_voltage_type()
        assert isinstance(volt, exago.OPFLOWGenBusVoltageType)
        assert volt == types[i]


@exago_test
def test_set_genbusvoltage_strings():
    '''Testing setting gen bus voltage type via strings'''
    opf = exago.OPFLOW()
    types = list(exago.OPFLOWGenBusVoltageType.__members__.items())
    # List of tuple pairs of stringy values and corresponding enum value
    print(types)
    for i in range(0, len(types)):
        opf.set_gen_bus_voltage_type(types[i][0])
        volt = opf.get_gen_bus_voltage_type()
        assert isinstance(volt, exago.OPFLOWGenBusVoltageType)
        assert volt == types[i][1]


@exago_test
def test_set_has_gensetpoint():
    '''Testing setting has gen set point'''
    opf = exago.OPFLOW()
    opf.set_has_gen_set_point(True)
    b = opf.get_has_gen_set_point()
    assert isinstance(b, bool)
    assert b is True
    opf.set_has_gen_set_point(False)
    b = opf.get_has_gen_set_point()
    assert isinstance(b, bool)
    assert b is False


@exago_test
def test_obj_func():
    '''Testing objective function value'''
    opf = exago.OPFLOW()
    path = exago.prefix()
    opf.read_mat_power_data(os.path.join(
        path, 'share', 'exago', 'datafiles', 'case9', 'case9mod.m'))
    opf.solve()
    obj = opf.get_objective()
    assert isinstance(obj, float)


@exago_test
def test_convergence_status():
    '''Testing convergence status'''
    opf = exago.OPFLOW()
    path = exago.prefix()
    opf.read_mat_power_data(os.path.join(
        path, 'share', 'exago', 'datafiles', 'case9', 'case9mod.m'))
    opf.solve()
    s = opf.get_convergence_status()
    assert isinstance(s, bool)
    assert s is True


@exago_test
def test_set_hiop_compute_mode():
    '''Testing hiop compute mode'''
    opf = exago.OPFLOW()
    opf.set_solver("HIOP")
    opf.set_model("POWER_BALANCE_HIOP")
    opf.set_hiop_compute_mode("CPU")
    mode = opf.get_hiop_compute_mode()
    assert mode == "CPU"


@exago_test
def test_set_hiop_verbosity():
    '''Testing setting hiop verbosity level'''
    opf = exago.OPFLOW()
    v = 3
    opf.set_hiop_verbosity_level(v)
    level = opf.get_hiop_verbosity_level()
    assert isinstance(level, int)
    assert level == v


@exago_test
def test_set_model_and_solver():
    '''Testing setting model'''
    opf = exago.OPFLOW()
    opf.set_solver("HIOP")
    solver = opf.get_solver()
    assert solver == "HIOP"
    opf.set_model("POWER_BALANCE_HIOP")
    model = opf.get_model()
    assert model == "POWER_BALANCE_HIOP"
    opf.set_hiop_verbosity_level(3)
    opf.set_hiop_compute_mode('cpu')
    path = exago.prefix()
    opf.read_mat_power_data(os.path.join(
        path, 'share', 'exago', 'datafiles', 'case9', 'case9mod.m'))
    opf.solve()
    obj = opf.get_objective()
    assert isinstance(obj, float)


@exago_test
def test_tolerance_opflow():
    '''Testing setting tolerance in opflow object'''
    opf = exago.OPFLOW()
    opf.set_tolerance(1e-6)
    tol = opf.get_tolerance()
    assert isinstance(tol, float)
    assert tol == 1e-6
    opf.set_tolerance(1e-5)
    tol = opf.get_tolerance()
    assert isinstance(tol, float)
    assert tol == 1e-5


@exago_test
def test_set_loadloss_penalty_opflow():
    '''Testing setting load loss penalty in opflow object'''
    opf = exago.OPFLOW()
    opf.set_loadloss_penalty(999)
    penalty = opf.get_loadloss_penalty()
    assert isinstance(penalty, float)
    assert penalty == 999


@exago_test
def test_ignore_lineflow_constraints():
    '''Testing toggling use lineflow constraints'''
    opf = exago.OPFLOW()
    opf.set_ignore_lineflow_constraints(True)
    b = opf.get_ignore_lineflow_constraints()
    assert isinstance(b, bool)
    assert b is True
    opf.set_ignore_lineflow_constraints(False)
    b = opf.get_ignore_lineflow_constraints()
    assert isinstance(b, bool)
    assert b is False


@exago_test
def test_include_loadloss():
    '''Testing toggling use loadloss'''
    opf = exago.OPFLOW()
    opf.set_has_loadloss(True)
    b = opf.get_has_loadloss()
    assert isinstance(b, bool)
    assert b is True
    opf.set_has_loadloss(False)
    b = opf.get_has_loadloss()
    assert isinstance(b, bool)
    assert b is False


@exago_test
def test_include_powerimbalance():
    '''Testing toggling use powerimbalance'''
    opf = exago.OPFLOW()
    opf.set_has_bus_power_imbalance(True)
    b = opf.get_has_bus_power_imbalance()
    assert isinstance(b, bool)
    assert b is True
    opf.set_has_bus_power_imbalance(False)
    b = opf.get_has_bus_power_imbalance()
    assert isinstance(b, bool)
    assert b is False


@exago_test
def test_use_agc():
    '''Testing toggling AGC'''
    opf = exago.OPFLOW()
    opf.set_use_agc(True)
    b = opf.get_use_agc()
    assert isinstance(b, bool)
    assert b is True
    opf.set_use_agc(False)
    b = opf.get_use_agc()
    assert isinstance(b, bool)
    assert b is False


@exago_test
def test_set_bus_powerimbalance_penalty_opflow():
    '''Testing setting power imbalance penalty in opflow object'''
    opf = exago.OPFLOW()
    opf.set_bus_power_imbalance_penalty(999)
    p = opf.get_bus_power_imbalance_penalty()
    assert isinstance(p, float)
    assert p == 999


@exago_test
def test_solve_opflow():
    '''Testing opflow solve'''
    opf = exago.OPFLOW()
    path = exago.prefix()
    opf.read_mat_power_data(os.path.join(
        path, 'share', 'exago', 'datafiles', 'case9', 'case9mod.m'))
    opf.solve()


@exago_test
def test_ps():
    '''Testing ps set up'''
    opf = exago.OPFLOW()
    path = exago.prefix()
    opf.read_mat_power_data(os.path.join(
        path, 'share', 'exago', 'datafiles', 'case9', 'case9mod.m'))
    opf.set_up_ps()


@exago_test
def test_ps_set_gen_power_limits():
    '''Testing ps set gen power limits'''
    opf = exago.OPFLOW()
    opf.set_solver("HIOP")
    opf.set_model("POWER_BALANCE_HIOP")
    opf.set_hiop_compute_mode('cpu')
    path = exago.prefix()
    opf.read_mat_power_data(os.path.join(
        path, 'share', 'exago', 'datafiles', 'case9', 'case9mod.m'))
    opf.set_up_ps()
    opf.ps_set_gen_power_limits(3, "1 ", 65, 0, exago_ignore, exago_ignore)
    opf.solve()


@exago_test
def test_solution_to_ps():
    '''Testing ps solution'''
    opf = exago.OPFLOW()
    opf.set_solver("HIOP")
    opf.set_model("POWER_BALANCE_HIOP")
    opf.set_hiop_compute_mode('cpu')
    path = exago.prefix()
    opf.read_mat_power_data(os.path.join(
        path, 'share', 'exago', 'datafiles', 'case9', 'case9mod.m'))
    opf.set_up_ps()
    opf.ps_set_gen_power_limits(3, "1 ", 65, 0, exago_ignore, exago_ignore)
    opf.solve()
    opf.solution_to_ps()


@exago_test
def test_get_gen_dispatch():
    '''Testing ps get gen dispatch'''
    opf = exago.OPFLOW()
    opf.set_solver("HIOP")
    opf.set_model("POWER_BALANCE_HIOP")
    opf.set_hiop_compute_mode('cpu')
    opf.set_hiop_verbosity_level(3)
    path = exago.prefix()
    opf.read_mat_power_data(os.path.join(
        path, 'share', 'exago', 'datafiles', 'case9', 'case9mod.m'))
    opf.set_up_ps()
    opf.ps_set_gen_power_limits(3, "1 ", 65, 0, exago_ignore, exago_ignore)
    opf.solve()
    opf.solution_to_ps()
    pg1set, qg1set = opf.get_gen_dispatch(1, "1 ")
    assert isinstance(pg1set, float)
    assert isinstance(qg1set, float)


@exago_test
def test_save_solution():
    '''Testing opflow saving solve to file'''
    opf = exago.OPFLOW()
    opf.set_solver("HIOP")
    opf.set_model("POWER_BALANCE_HIOP")
    opf.set_hiop_compute_mode('cpu')
    path = exago.prefix()
    opf.read_mat_power_data(os.path.join(
        path, 'share', 'exago', 'datafiles', 'case9', 'case9mod.m'))
    opf.solve()
    filepath = os.path.join(path, '..', 'tmp.output')
    opf.save_solution(exago.CSV, filepath)


@exago_test
def test_finalize():
    '''Test calling ExaGOFinalize'''
    exago.finalize()


exago_run_all_tests(__file__)
