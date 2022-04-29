from check_preconditions import check_preconditions
from test_utilities import *
import mpi4py.rc
mpi4py.rc.threads = False
from mpi4py import MPI  # noqa
import exago  # noqa
check_preconditions()


@exago_test
def test_initialize():
    '''Test call ExaGOInitialize'''
    comm = MPI.COMM_WORLD
    exago.initialize("pflow", comm)


@exago_test
def test_get_datafile_path():
    '''Test retrieving datafile path'''
    path = exago.prefix()
    assert isinstance(path, str)


@exago_test
def test_creating_pflow():
    '''Testing creation of pflow object'''
    pflow = exago.PFLOW()
    print(pflow)


@exago_test
def test_solve():
    '''Testing pflow solve'''
    pflow = exago.PFLOW()
    path = exago.prefix()
    pflow.read_mat_power_data(os.path.join(
        path, 'share', 'exago', 'datafiles', 'case9', 'case9mod.m'))
    pflow.solve()


@exago_test
def test_finalize():
    '''Test call ExaGOFinalize'''
    exago.finalize()


exago_run_all_tests(__file__)