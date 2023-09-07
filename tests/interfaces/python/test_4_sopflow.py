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

# SCOPFLOW has serious destructor problems, so it can't do this
# @pytest.mark.nocomm
# @pytest.mark.MPI
# def test_creating_scopflow():
#     '''Testing creation of sopflow object'''
#    opf = exago.SOPFLOW()

# -------------------------------------------------------------
# run_multicontingency
# -------------------------------------------------------------


def run_multicontingency(solver, nscen, ncont):
    sopf = exago.SOPFLOW()
    path = exago.prefix()

    sopf.set_tolerance(1.0E-03)

    sopf.set_network_data(os.path.join(
        path, 'share', 'exago', 'datafiles', 'case9', 'case9mod_gen3_wind.m'))

    sopf.set_num_scenarios(nscen)
    sopf.set_scenario_data(os.path.join(
        path, 'share', 'exago', 'datafiles', 'case9', '10_scenarios_9bus.csv'),
        exago.ScenarioFileInputFormat.NATIVE_SINGLEPERIOD,
        exago.ScenarioUncertaintyType.WIND
    )

    sopf.enable_multi_contingency(True)
    sopf.set_contingency_data(os.path.join(
        path, 'share', 'exago', 'datafiles', 'case9', 'case9.cont'),
        exago.ContingencyFileInputFormat.NATIVE
    )
    sopf.set_num_contingencies(ncont)

    sopf.set_solver('IPOPT')

    sopf.setup()
    sopf.solve()

    assert sopf.get_convergence_status()

    obj = sopf.get_base_objective()
    assert isinstance(obj, float)

    obj = sopf.get_total_objective()
    assert isinstance(obj, float)

    n = sopf.get_num_iterations()
    assert isinstance(n, int)

    oname = os.path.join(tempfile.gettempdir(), "sopflow_test_solution.csv")
    sopf.save_solution(0, exago.OutputFormat.CSV, oname)
    assert os.path.exists(oname)
    os.unlink(oname)

    oname = os.path.join(tempfile.gettempdir(), "sopflow_test_solution.m")
    sopf.save_solution(0, exago.OutputFormat.MATPOWER, oname)
    assert os.path.exists(oname)
    os.unlink(oname)

    oname = os.path.join(tempfile.gettempdir(), "sopflow_test_solution_dir")
    sopf.save_solution_all(exago.OutputFormat.MATPOWER, oname)
    assert os.path.exists(oname)
    shutil.rmtree(oname)


# From sopflow_multicontingency.toml:
# description =
# 'datafiles/case9/case9mod_gen3_wind.m IPOPT 3scen contingencies SOPFLOW'
@pytest.mark.nocomm
@pytest.mark.MPI
def test_multicontingency_1():
    run_multicontingency('HIOP', 3, 1)


@pytest.mark.nocomm
@pytest.mark.MPI
def test_multicontingency_2():
    run_multicontingency('HIOP', 3, 2)


@pytest.mark.nocomm
@pytest.mark.MPI
def test_multicontingency_3():
    run_multicontingency('HIOP', 3, 3)


@pytest.mark.nocomm
@pytest.mark.MPI
def test_multicontingency_4():
    run_multicontingency('EMPAR', 3, 1)


@pytest.mark.nocomm
@pytest.mark.MPI
def test_multicontingency_5():
    run_multicontingency('EMPAR', 3, 2)


@pytest.mark.nocomm
@pytest.mark.MPI
def test_multicontingency_6():
    run_multicontingency('EMPAR', 3, 3)


@pytest.mark.nocomm
@pytest.mark.MPI
def test_multicontingency_7():
    run_multicontingency('HIOP', 3, 6)


@pytest.mark.nocomm
@pytest.mark.MPI
def test_multicontingency_8():
    sopf = exago.SOPFLOW()
    path = exago.prefix()

    nscen = 3
    ncont = 6

    sopf.set_tolerance(1.0E-03)

    sopf.set_network_data(os.path.join(
        path, 'share', 'exago', 'datafiles', 'case9', 'case9mod_gen3_wind.m'))

    sopf.set_num_scenarios(nscen)
    sopf.set_scenario_data(os.path.join(
        path, 'share', 'exago', 'datafiles', 'case9', '10_scenarios_9bus.csv'),
        exago.ScenarioFileInputFormat.NATIVE_SINGLEPERIOD,
        exago.ScenarioUncertaintyType.WIND
    )

    sopf.set_initialization_type(
        exago.OPFLOWInitializationType.OPFLOWINIT_ACPF)
    sopf.set_gen_bus_voltage_type('VARIABLE_WITHIN_BOUNDS')

    sopf.enable_multi_contingency(True)
    sopf.set_contingency_data(os.path.join(
        path, 'share', 'exago', 'datafiles', 'case9', 'case9.cont'),
        exago.ContingencyFileInputFormat.NATIVE
    )
    sopf.set_num_contingencies(ncont)

    sopf.set_solver('IPOPT')

    sopf.setup()
    sopf.solve()

    assert sopf.get_convergence_status()

    obj = sopf.get_base_objective()
    assert isinstance(obj, float)

    obj = sopf.get_total_objective()
    assert isinstance(obj, float)

    n = sopf.get_num_iterations()
