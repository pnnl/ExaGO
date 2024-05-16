# test only run once things here
import pytest
from check_preconditions import check_preconditions
check_preconditions()
import exago

@pytest.mark.nocomm
@pytest.mark.MPI
def test_finalize():
    '''Test calling ExaGOFinalize'''
    exago.finalize()