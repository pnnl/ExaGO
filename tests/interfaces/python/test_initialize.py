from check_preconditions import check_preconditions
from test_utilities import *
import mpi4py.rc
mpi4py.rc.threads = False
import exago  # noqa
check_preconditions()


@exago_test
def test_initialize():
    '''Test call ExaGOInitialize without a communicator'''
    exago.initialize("no-comm")


@exago_test
def test_finalize():
    '''Test call ExaGOFinalize'''
    exago.finalize()


exago_run_all_tests(__file__)
