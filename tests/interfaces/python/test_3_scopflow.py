import os
import pytest
from check_preconditions import check_preconditions
import mpi4py.rc
mpi4py.rc.threads = False
from mpi4py import MPI  # noqa
import exago  # noqa

check_preconditions()

exago_ignore = -1000000


@pytest.mark.nocomm
@pytest.mark.MPI
def test_creating_scopflow():
    '''Testing creation of scopflow object'''
    opf = exago.SCOPFLOW()


@pytest.mark.nocomm
@pytest.mark.MPI
def test_set_initialization_ints():
    '''Testing setting initialization type via integers'''
    scopf = exago.SCOPFLOW()
    opf = exago.OPFLOW()
    types = opf.get_initialization_types()
    del opf
    # List of possible enum values (object list)
    print(types)
    for i in range(0, len(types)):
        scopf.set_initialization_type(types[i])


@pytest.mark.nocomm
@pytest.mark.MPI
def test_set_initialization_strings():
    '''Testing setting initialization type via strings'''
    scopf = exago.SCOPFLOW()
    types = list(exago.OPFLOWInitializationType.__members__.items())
    # List of tuple pairs of stringy values and corresponding enum value
    print(types)
    for i in range(0, len(types)):
        scopf.set_initialization_type(types[i][0])


@pytest.mark.nocomm
@pytest.mark.MPI
def test_set_genbusvoltage_ints():
    '''Testing setting gen bus voltage type via ints'''
    scopf = exago.SCOPFLOW()
    opf = exago.OPFLOW()
    types = opf.get_gen_bus_voltage_types()
    del opf
    # List of possible enum values (object list)
    print(types)
    for i in range(0, len(types)):
        scopf.set_gen_bus_voltage_type(types[i])


@pytest.mark.nocomm
@pytest.mark.MPI
def test_set_genbusvoltage_strings():
    '''Testing setting gen bus voltage type via strings'''
    scopf = exago.SCOPFLOW()
    types = list(exago.OPFLOWGenBusVoltageType.__members__.items())
    # List of tuple pairs of stringy values and corresponding enum value
    print(types)
    for i in range(0, len(types)):
        scopf.set_gen_bus_voltage_type(types[i][0])


@pytest.mark.nocomm
@pytest.mark.MPI
def test_tolerance_scopflow():
    '''Testing setting tolerance in scopflow object'''
    scopf = exago.SCOPFLOW()
    scopf.set_tolerance(1e-6)
    tol = scopf.get_tolerance()
    assert isinstance(tol, float)
    assert tol == 1e-6
    scopf.set_tolerance(1e-5)
    tol = scopf.get_tolerance()
    assert isinstance(tol, float)
    assert tol == 1e-5


@pytest.mark.nocomm
@pytest.mark.MPI
def test_set_network():
    '''Test reading network data'''
    scopf = exago.SCOPFLOW()
    path = exago.prefix()
    scopf.set_network_data(
        os.path.join(path, 'share', 'exago',
                     'datafiles', 'case9', 'case9mod.m'))


@pytest.mark.nocomm
@pytest.mark.MPI
def test_set_contingency_data():
    '''Test reading contingency data'''
    scopf = exago.SCOPFLOW()
    path = exago.prefix()
    scopf.set_network_data(
        os.path.join(path, 'share', 'exago',
                     'datafiles', 'case9', 'case9mod.m'))
    scopf.set_contingency_data(
        os.path.join(path, 'share', 'exago',
                     'datafiles', 'case9', 'case9.cont'),
        exago.ContingencyFileInputFormat.NATIVE
    )


@pytest.mark.nocomm
@pytest.mark.MPI
def test_complete_scopflow_solve():
    '''Testing scopflow solve'''
    scopf = exago.SCOPFLOW()
    path = exago.prefix()
    scopf.set_network_data(os.path.join(
        path, 'share', 'exago', 'datafiles', 'case9', 'case9mod.m'))
    scopf.set_contingency_data(os.path.join(
        path, 'share', 'exago', 'datafiles', 'case9', 'case9.cont'),
        exago.ContingencyFileInputFormat.NATIVE
    )
    scopf.set_tolerance(1.0E-03)
    scopf.set_solver('IPOPT')
    scopf.set_model('GENRAMP')
    scopf.set_initialization_type(
        exago.OPFLOWInitializationType.OPFLOWINIT_ACPF
    )
    scopf.set_num_contingencies(-1)
    scopf.set_subproblem_model('POWER_BALANCE_POLAR')
    scopf.set_subproblem_solver('IPOPT')
    scopf.set_gen_bus_voltage_type('VARIABLE_WITHIN_BOUNDS')
    scopf.enable_multi_period(False)
    scopf.enable_power_imbalance_variables(True)
    scopf.set_verbosity_level(3)
    scopf.set_mode(0)
    scopf.set_up()
    scopf.solve()
