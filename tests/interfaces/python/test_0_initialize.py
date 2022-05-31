# test initializing exago only once
import os
import pytest
from check_preconditions import check_preconditions
import mpi4py.rc
mpi4py.rc.threads = False
from mpi4py import MPI  # noqa
import exago  # noqa


check_preconditions()


@pytest.mark.nocomm
def test_initialize_nocomm():
    '''Test call ExaGOInitialize without a communicator'''
    exago.initialize("no-comm")


@pytest.mark.MPI
def test_initialize_MPI():
    '''Test call ExaGOInitialize'''
    comm = MPI.COMM_WORLD
    exago.initialize("test", comm)
