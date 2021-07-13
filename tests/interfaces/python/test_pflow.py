from exago.pflow import PFLOW
from exago import config
from test_utilities import *
from ctypes import *
import logging
from check_preconditions import check_preconditions
check_preconditions()
   

@exago_test
def test_pflow_create():
    '''Test creation of powerflow object'''
    # Since __del__ calls MPI Finalize, we only want this called once
    # We overrite in each test to avoid calling this until the end.
    # In a non-testing environment this shouldn't be an issue(?)
    pf = PFLOW()
    pf.dont_finalize()


@exago_test
def test_pflow_read_mat_power_data():
    '''Test powerflow reading data in mat power file format'''
    pf = PFLOW()
    pf.dont_finalize()
    pf.read_mat_power_data(
        os.path.join(config.prefix(), 'share', 'exago', 'datafiles', 'case9', 'case9mod.m'))


@exago_test
def test_pflow_solve_case9():
    '''Test solving a powerflow problem'''
    pf = PFLOW()
    pf.read_mat_power_data(
        os.path.join(config.prefix(), 'share', 'exago', 'datafiles', 'case9', 'case9mod.m'))
    pf.solve()

exago_run_all_tests(__file__)
