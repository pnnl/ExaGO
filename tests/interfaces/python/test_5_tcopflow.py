import tempfile
import os
import shutil
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
def test_creating_tcopflow():
    '''Testing creation of tcopflow object'''
    tcopf = exago.TCOPFLOW()


@pytest.mark.nocomm
@pytest.mark.MPI
def test_get_prefix():
    '''Test retrieving datafile path'''
    path = exago.prefix()
    assert isinstance(path, str)


''' example
@pytest.mark.nocomm
@pytest.mark.MPI
def test_set_model_and_solver():
    <put three single quotes here> Testing setting model <put three single quotes here>
    opf = exago.OPFLOW()
    opf.set_solver("HIOP")
    solver = opf.get_solver()
    assert solver == "HIOP"
    opf.set_model("POWER_BALANCE_HIOP")
    model = opf.get_model()
    assert model == "POWER_BALANCE_HIOP"
    path = exago.prefix()
    opf.read_mat_power_data(os.path.join(
        path, 'share', 'exago', 'datafiles', 'case9', 'case9mod.m'))
'''


@pytest.mark.nocomm
@pytest.mark.MPI
def test_set_tolerance():
    '''test get/set tolerance for tcopflow'''
    tcopf = exago.TCOPFLOW()
    tcopf.set_tolerance(-6)
    tolerance = tcopf.get_tolerance()
    assert isinstance(tolerance, float)
    assert tolerance == -6
