# test only run once things here
import exago
import pytest
from check_preconditions import check_preconditions
check_preconditions()


@pytest.mark.nocomm
@pytest.mark.MPI
def test_finalize():
    '''Test calling ExaGOFinalize'''
    exago.finalize()
