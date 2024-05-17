from check_preconditions import check_preconditions
import os
import pytest
import mpi4py.rc
mpi4py.rc.threads = False
from mpi4py import MPI  # noqa
check_preconditions()
import exago  # noqa


@pytest.mark.nocomm
@pytest.mark.MPI
def test_get_datafile_path():
    '''Test retrieving datafile path'''
    path = exago.prefix()
    assert isinstance(path, str)


@pytest.mark.nocomm
@pytest.mark.MPI
def test_creating_pflow():
    '''Testing creation of pflow object'''
    pflow = exago.PFLOW()
    print(pflow)


@pytest.mark.nocomm
@pytest.mark.MPI
def test_solve():
    '''Testing pflow solve'''
    pflow = exago.PFLOW()
    path = exago.prefix()
    pflow.read_mat_power_data(os.path.join(
        path, 'share', 'exago', 'datafiles', 'case9', 'case9mod.m'))
    pflow.solve()
