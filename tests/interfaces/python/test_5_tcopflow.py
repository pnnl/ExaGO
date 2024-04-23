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

''' cannot run just creating tcopflow object bc opflow sub objects
   won't be allocated (done in TCOPFLOWSetUp) so destructor seg faults
def test_creating_tcopflow():
    tcopf = exago.TCOPFLOW()
'''
def run_tcopflow(solver):
    tcopf = exago.TCOPFLOW()
    path = exago.prefix()

    tcopf.set_tolerance(1.0E-03)
    tcopf.set_network_data(os.path.join(
        path, 'share', 'exago', 'datafiles', 'case9', 'case9mod_gen3_wind.m'))
    tcopf.set_solver(solver)
    #set_time_step_and_duration
    #set_load_profiles
    #set_wind_gen_profiles

    tcopf.setup()
    tcopf.solve()

    assert tcopf.get_convergence_status()

    obj = tcopf.get_objective()
    assert isinstance(obj, float)

    n = tcopf.get_num_iterations()
    assert isinstance(n, int)

    #get_solution
    #print_solution

    oname = os.path.join(tempfile.gettempdir(), "tcopflow_test_solution.csv")
    tcopf.save_solution(0, exago.OutputFormat.CSV, oname)
    assert os.path.exists(oname)
    os.unlink(oname)

    oname = os.path.join(tempfile.gettempdir(), "tcopflow_test_solution.m")
    tcopf.save_solution(0, exago.OutputFormat.MATPOWER, oname)
    assert os.path.exists(oname)
    os.unlink(oname)

    oname = os.path.join(tempfile.gettempdir(), "tcopflow_test_solution_dir")
    tcopf.save_solution_all(exago.OutputFormat.MATPOWER, oname)
    assert os.path.exists(oname)
    shutil.rmtree(oname)


# From sopflow_multicontingency.toml:
# description =
# 'datafiles/case9/case9mod_gen3_wind.m IPOPT 3scen contingencies SOPFLOW'
@pytest.mark.nocomm
@pytest.mark.MPI
def test_tcopflow_1():
    run_tcopflow('IPOPT')



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
