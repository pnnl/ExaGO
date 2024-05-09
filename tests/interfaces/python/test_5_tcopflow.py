from mpi4py import MPI
import tempfile
import os
import shutil
import pytest
import mpi4py.rc
mpi4py.rc.threads = False
from check_preconditions import check_preconditions
check_preconditions()
import exago

check_preconditions()


def run_tcopflow(solver):
    tcopf = exago.TCOPFLOW()
    path = exago.prefix()

    tcopf.set_tolerance(1.0E-03)
    tcopf.set_network_data(os.path.join(
        path, 'share', 'exago', 'datafiles', 'case9', 'case9mod_gen3_wind.m'))
    tcopf.set_solver(solver)

    tcopf.setup()
    tcopf.solve()

    assert tcopf.get_convergence_status()

    obj = tcopf.get_objective()
    assert isinstance(obj, float)

    n = tcopf.get_num_iterations()
    assert isinstance(n, int)


@pytest.mark.nocomm
@pytest.mark.MPI
def test_tcopflow_1():
    run_tcopflow('IPOPT')
