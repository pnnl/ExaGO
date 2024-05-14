## ExaGO Python Bindings

### Overview

The ExaGO Python bindings use an object-oriented API slightly different from the
C++ API.
The C++ API uses the application type in uppercase as the prefix for its
methods, where they are native methods in Python.
To view the implementation of this wrapper, see [/interfaces/python/exago_python.cpp](../interfaces/python/exago_python.cpp).

For example, solving an OPF from C++ might look like this:
```cpp
#include <opflow.h>
#include <exago_config.h>

static char help[] = "User example calling OPFLOW.\n";
static char appname[] = "opflow";

int main(int argc, char** argv) {

  /* Initialize ExaGO application */
  PetscErrorCode ierr;
  OPFLOW opflow;
  MPI_Comm comm = MPI_COMM_WORLD;
  ierr = ExaGOInitialize(comm, &argc, &argv, appname, help);

  /* Create OPFLOW object */
  ierr = OPFLOWCreate(comm, &opflow);
  ExaGOCheckError(ierr);

  /* Read network data */
  ierr = OPFLOWReadMatPowerData(opflow, "datafiles/case9/case9mod.m");
  ExaGOCheckError(ierr);

  /* Solve */
  ierr = OPFLOWSolve(opflow);
  ExaGOCheckError(ierr);

  /* Print solution */
  ierr = OPFLOWPrintSolution(opflow);
  ExaGOCheckError(ierr);

  /* Destroy OPFLOW object */
  ierr = OPFLOWDestroy(&opflow);
  ExaGOCheckError(ierr);

  /* Clean up resources */
  ExaGOFinalize();

  return 0;
}
```

While a version of the same driver that uses the Python bindings looks like
this:

```python
import exago
exago.initialize("opflow")
opf = exago.OPFLOW()
opf.read_mat_power_data('datafiles/case9/case9mod.m')
opf.solve()
opf.print_solution()
del opf
exago.finalize()
```

Note that the Python bindings use snake case in accordance with standard Python
naming conventions.
Translating between the C++ API and the Python API should be relatively
straightforward.  
ExaGO Python instances must be destroyed, using `del`, before
calling `exago.finalize()`, like calling `*Destroy()` with the C++
API.  Failure to do so will cause segmentation faults or other memory
errors.  

***If you identify components of the C++ API that you need to call from Python,
please [open an issue on our issues page](https://gitlab.pnnl.gov/exasgd/frameworks/exago/-/issues).***

### Building with MPI

ExaGO depends on mpi4py when running through the python interface.

When running exago, you may need to disable threading through mpi4py before importing exago:

```python
import mpi4py.rc
mpi4py.rc.threads = False
from mpi4py import MPI
import exago
comm = MPI.COMM_WORLD
exago.initialize("app", comm)
# ...
exago.finalzie()
```

Additionally linting tools may re-order your imports, and so you may need to add appropriate comments in order to avoid this:

```python
import mpi4py.rc
mpi4py.rc.threads = False
from mpi4py import MPI # noqa
import exago # noqa
```

### Bindings Tables

#### ExaGO

| C++ API | Python API | Notes |
|---|---|---|
| `ExaGOInitialize` | `exago.initialize` |  |
| `ExaGOFinalize` | `exago.finalize` |  |
|  | `exago.prefix` |  Returns the path to the installation directory of ExaGO (for finding mat power data files, etc) |
| `OutputFormat` enum | `exago.OutputFormat` enum | Output format type for functions like `save_solution`. Possible values for this enum can be found in the [next section](#enums). |

#### PFLOW

| C++ API | Python API | Notes |
|---|---|---|
| `PFLOW` | `exago.PFLOW` class | From here below, `pflow` objects are instances of the Python `exago.PFLOW` class. The same is the case for opflow, scopflow, tcopflow, and sopflow. |
| `PFLOWReadMatPowerData` | `pflow.read_mat_power_data` |  |
| `PFLOWSolve` | `pflow.solve` |  |

#### OPFLOW

The following table assumes `opflow = exago.OPFLOW()`.

| C++ API | Python API | Notes |
|---|---|---|
| `OPFLOW` | `exago.OPFLOW` class |  |
| `OPFLOWObjectiveType` enum | `exago.OPFLOWObjectiveType` enum | More details and possible values for this enum can be found in the [next section](#enums). |
| `OPFLOWInitializationType` enum | `exago.OPFLOWInitializationType` enum | More details and possible values for this enum can be found in the [next section](#enums). |
| `OPFLOWGenBusVoltageType` enum | `exago.OPFLOWGenBusVoltageType` enum | More details and possible values for this enum can be found in the [next section](#enums). |
| `OPFLOWSetObjectiveType` | `opflow.set_objective_type` | |
| `OPFLOWSetInitializationType` | `opflow.set_initialization_type` | |
| `OPFLOWSetGenBusVoltageType` | `opflow.set_gen_bus_voltage_type` | |
| `OPFLOWSetModel` | `opflow.set_model` | options are "PBPOLRAJAHIOP", "POWER_BALANCE_HIOP", and "POWER_BALANCE_POLAR" |
| `OPFLOWSetSolver` | `opflow.set_solver` | options are "IPOPT", "HIOP", and "HIOPSPARSE"  |
| `OPFLOWHasGenSetPoint` | `opflow.set_has_gen_set_point` | |
| `OPFLOWSetHIOPComputeMode` | `opflow.set_hiop_compute_mode` | options are "CPU" or "GPU"|
| `OPFLOWSetHIOPMemSpace` | `opflow.set_hiop_mem_space` | options are "DEFAULT", "HOST", "UM", and "DEVICE" |
| `OPFLOWHasLoadLoss` | `opflow.set_has_loadloss` | |
| `OPFLOWIgnoreLineflowConstraints` | `opflow.set_ignore_lineflow_constraints` | |
| `OPFLOWHasBusPowerImbalance` | `opflow.set_has_bus_power_imbalance` | |
| `OPFLOWUseAGC` | `opflow.set_use_agc` | |
| `OPFLOWSetHIOPVerbosityLevel` | `opflow.set_hiop_verbosity_level` | integer between 0 and 10 |
| `OPFLOWSetLoadLossPenalty` | `opflow.set_loadloss_penalty` | |
| `OPFLOWSetBusPowerImbalancePenalty` | `opflow.set_bus_power_imbalance_penalty` | |
| `OPFLOWSetTolerance` | `opflow.set_tolerance` | |
| `OPFLOWSetWeight` | `opflow.set_weight` | |
| `PSSetGenPowerLimits` | `opflow.ps_set_gen_power_limits` | |
| `OPFLOWGetTolerance` | `opflow.get_tolerance` | |
| `OPFLOWGetHIOPComputeMode` | `opflow.get_hiop_compute_mode` | |
| `OPFLOWGetHIOPMemSpace` | `opflow.get_hiop_mem_space` | |
| `OPFLOWGetModel` | `opflow.get_model` | |
| `OPFLOWGetSolver` | `opflow.get_solver` | |
| `OPFLOWGetConvergenceStatus` | `opflow.get_convergence_status` | |
| `OPFLOWGetObjectiveType` | `opflow.get_objective_type` | |
| `OPFLOWGetInitializationType` | `opflow.get_initialization_type` | |
| `OPFLOWGetGenBusVoltageType` | `opflow.get_gen_bus_voltage_type` | |
| `OPFLOWGetHasGenSetPoint` | `opflow.get_has_gen_set_point` | |
| `OPFLOWGetLoadlossPenalty` | `opflow.get_loadloss_penalty` | |
| `OPFLOWGetIgnoreLineflowConstraints` | `opflow.get_ignore_lineflow_constraints` | |
| `OPFLOWGetHasLoadloss` | `opflow.get_has_loadloss` | |
| `OPFLOWGetHasBusPowerImbalance` | `opflow.get_has_bus_power_imbalance` | |
| `OPFLOWGetUseAGC` | `opflow.get_use_agc` | |
| `OPFLOWGetHIOPVerbosityLevel` | `opflow.get_hiop_verbosity_level` | |
| `OPFLOWGetBusPowerImbalancePenalty` | `opflow.get_bus_power_imbalance_penalty` | |
| `PSGetGenDispatch` | `opflow.get_gen_dispatch` | |
| `OPFLOWGetObjectiveTypes` | `opflow.get_objective_types` | |
| `OPFLOWGetInitializationTypes` | `opflow.get_initialization_types` | |
| `OPFLOWGetGenBusVoltageTypes` | `opflow.get_gen_bus_voltage_types` | |
| `OPFLOWGetObjective` | `opflow.get_objective` | |
| `OPFLOWSolve` | `opflow.solve` | |
| `OPFLOWPrintSolution` | `opflow.print_solution` | |
| `OPFLOWSaveSolution` | `opflow.save_solution` | |
| `OPFLOWReadMatPowerData` | `opflow.read_mat_power_data` | |
| `OPFLOWSolutionToPS` | `opflow.solution_to_ps` | |
| `OPFLOWSetUpPS` | `opflow.set_up_ps` | |
| `OPFLOWSkipOptions` | `opflow.skip_options | |
| `OPFLOWSetLinesMonitored` | | Implemented as two different methods: | |
| |`opflow.set_lines_monitored([...])` | Specify a list of line kvlevels (type `float`) to monitor |
| |`opflow.set_lines_monitored(n, "file")` | Read `n` line kvlevels from a file (`n=-1` for all). |

#### SCOPFLOW

The following table assumes `scopflow = exago.SCOPFLOW()`.

| C++ API | Python API | Notes |
|---|---|---|
| `SCOPFLOW` | `exago.SCOPFLOW` class |  |
| `ContingencyFileInputFormat` enum | `exago.ContingencyFileInputFormat` enum | More details and possible values for this enum can be found in the [next section](#enums).  Currently, this is only be used as a Python enum. A string representation is not available. |
| `SCOPFLOWSetModel` | `set_model` |  |
| `SCOPFLOWSetNetworkData` | `set_network_data` |  |
| `SCOPFLOWSetLoadProfiles` | `set_load_profiles` | |
| `SCOPFLOWSetNumContingencies` | `set_num_contingencies` |  |
| `SCOPFLOWSetContingencyData` | `set_contingency_data` |  |
| `SCOPFLOWSetPLoadData` | `set_pload_data` |  |
| `SCOPFLOWSetQLoadData` | `set_qload_data` |  |
| `SCOPFLOWSetWindGenProfile` | `set_wind_gen_profile` |  |
| `SCOPFLOWSetTimeStep` | `set_time_step` |  |
| `SCOPFLOWSetDuration` | `set_duration` |  |
| `SCOPFLOWSetTimeStepandDuration` | `set_time_step_and_duration` |  |
| `SCOPFLOWSetTolerance` | `set_tolerance` |  |
| `SCOPFLOWSetVerbosityLevel` | `set_verbosity_level` |  |
| `SCOPFLOWSetComputeMode` | `set_compute_mode` |  |
| `SCOPFLOWSetSolver` | `set_solver` |  |
| `SCOPFLOWSetSubproblemModel` | `set_subproblem_model` |  |
| `SCOPFLOWSetSubproblemSolver` | `set_subproblem_solver` |  |
| `SCOPFLOWSetInitilizationType` | `set_initialization_type` |  |
| `SCOPFLOWSetGenBusVoltageType` | `set_gen_bus_voltage_type` |  |
| `SCOPFLOWEnableMultiPeriod` | `enable_multi_period` |  |
| `SCOPFLOWEnablePowerImbalanceVariables` | `enable_power_imbalance_variables` |  |
| `SCOPFLOWIgnoreLineflowConstraints` | `ignore_lineflow_constraints` |  |
| `SCOPFLOWGetTolerance` | `get_tolerance` |  |
| `SCOPFLOWGetNumIterations` | `get_num_iterations` |  |
| `SCOPFLOWGetConvergenceStatus` | `get_convergence_status` |  |
| `SCOPFLOWGetTotalObjective` | `get_total_objective` |  |
| `SCOPFLOWGetBaseObjective` | `get_base_objective` |  |
| `SCOPFLOWSetUp` | `set_up` |  |
| `SCOPFLOWSolve` | `solve` |  |
| `SCOPFLOWPrintSolution` | `print_solution` |  |
| `SCOPFLOWSaveSolution` | `save_solution` |  |
| `SCOPFLOWSaveSolutionDefault` | `save_solution_default` |  |
| `SCOPFLOWSaveSolutionAll` | `save_solution_all` |  |
| `SCOPFLOWSaveSolutionAllDefault` | `save_solution_all_default` |  |

#### SOPFLOW

The following table assumes `sopflow = exago.SOPFLOW()`.

| C++ API | Python API | Notes |
|---|---|---|
| `SCOPFLOW` | `exago.SCOPFLOW` class |  |
| `ScenarioFileInputFormat` enum | `ScenarioFileInputFormat` enum  | More details and possible values for this enum can be found in the [next section](#enums).  Currently, this is only be used as a Python enum. A string representation is not available. |
| `ScenarioUncertaintyType` enum | `ScenarioUncertaintyType` enum  | More details and possible values for this enum can be found in the [next section](#enums).  Currently, this is only be used as a Python enum. A string representation is not available. |
| ` SOPFLOWSetModel` | `set_model` |  |
| ` SOPFLOWSetNetworkData` | `set_network_data` |  |
| ` SOPFLOWSetContingencyData` | `set_contingency_data` |  |
| ` SOPFLOWSetNumContingencies` | `set_num_contingencies` |  |
| ` SOPFLOWSetScenarioData` | `set_scenario_data` |  |
| ` SOPFLOWSetNumScenarios` | `set_num_scenarios` |  |
| ` SOPFLOWSetWindGenProfile` | `set_wind_gen_profile` |  |
| ` SOPFLOWSetTimeStepandDuration` | `set_time_step_and_duration` |  |
| ` SOPFLOWSetTolerance` | `set_tolerance` |  |
| ` SOPFLOWSetSubproblemVerbosityLevel` | `set_subproblem_verbosity_level` |  |
| ` SOPFLOWSetSubproblemComputeMode` | `set_subproblem_compute_mode` |  |
| ` SOPFLOWSetSubproblemModel` | `set_subproblem_model` |  |
| ` SOPFLOWSetSubproblemSolver` | `set_subproblem_solver` |  |
| ` SOPFLOWSetSolver` | `set_solver` |  |
| ` SOPFLOWSetInitializationType` | `set_initialization_type ` |  |
| ` SOPFLOWSetGenBusVoltageType` | `set_gen_bus_voltage_type` |  |
| ` SOPFLOWSetLoadProfiles` | `set_load_profiles` |  |
| ` SOPFLOWSetLoadProfiles` | `set_ignore_lineflow_constraints` |  |
| ` SOPFLOWEnableMultiContingency` | `enable_multi_contingency` |  |
| ` SOPFLOWFlattenContingencies` | `flatten_contingencies` | |
| ` SOPFLOWGetNumScenarios` | `get_num_scenarios` |  |
| ` SOPFLOWGetNumIterations` | `get_num_iterations` |  |
| ` SOPFLOWGetConvergenceStatus` | `get_convergence_status` |  |
| ` SOPFLOWGetTotalObjective` | `get_total_objective` |  |
| ` SOPFLOWGetConvergenceStatus` | `get_converged_status` |  |
| ` SOPFLOWGetTolerance` | `get_tolerance` |  |
| ` SOPFLOWSetUp` | `setup` |  |
| ` SOPFLOWSolve` | `solve` |  |
| ` SOPFLOWPrintSolution` | `print_solution` | |
| ` SOPFLOWSaveSolution` | `save_solution` | |
| ` SOPFLOWSaveSolutionAll` | `save_solution_all` |  |

#### TCOPFLOW

The following table assumes `tcopflow = exago.TCOPFLOW()`.

| C++ API | Python API | Notes |
|---|---|---|
| `TCCOPFLOW` | `exago.TCCOPFLOW` class |  |
| ` TCOPFLOWSetNetworkData` | `set_network_data` |  |
| ` TCOPFLOWSetSolver` | `set_solver` |  |
| ` TCOPFLOWSetTolerance` | `set_tolerance` |  |
| ` TCOPFLOWSetup` | `setup` |  |
| ` TCOPFLOWSolve` | `solve` |  |
| ` TCOPFLOWGetConvergenceStatus` | `get_convergence_status` |  |
| ` TCOPFLOWGetObjective` | `get_objective` |  |
| ` TCOPFLOWGetNumIterations` | `get_num_iterations` |  |

### Enums

`OutputFormat` is the enum that specifies the output format in functions like `opflow.save_solution`.

Instances can be constructed directly through the `exago` library (i.e. `exago.OutputFormat.CSV`). 

OPFLOW has several type settings that are represented by enums: `OPFLOWObjectiveType`, `OPFLOWInitializationType`, and `OPFLOWGenBusVoltageType`. 

The possible values for the current enums are as follows:

OutputFormat
* CSV 
* MATPOWER

OPFLOWObjectiveType
* MIN_GEN_COST
* MIN_GENSETPOINT_DEVIATION
* NO_OB

OPFLOWInitializationType
* OPFLOWINIT_FROMFILE 
* OPFLOWINIT_MIDPOINT
* OPFLOWINIT_ACPF
* OPFLOW_FLATSTART

OPFLOWGenBusVoltageType
* VARIABLE_WITHIN_BOUNDS
* FIXED_WITHIN_QBOUNDS
* FIXED_AT_SETPOINT

ContingencyFileInputFormat
* NATIVE
* PSSE

ScenarioFileInputFormat
* NATIVE_SINGLEPERIOD
* NATIVE_MULTIPERIOD

ScenarioUncertaintyType
* NONE
* WIND
* LOAD

Note: getters for the possible values of ContingencyFileInputeFormat, ScenarioFileInputFormat, ScenarioUncertaintyType are not currently available

Instances can be constructed directly through the `exago` library (i.e. `exago.OPFLOWObjectiveType.MIN_GEN_COST` or `exago.MIN_GEN_COST`)
Setter functions for these OPFLOW Type configurations can take an integer, an instance of the exago enum, or a string that describes the enum (i.e. 'MIN_GEN_COST'). The rest currently only accept an instance of the exago enum.

The possible values for these OPFLOW Type enums can be retrieved through `opflow.get_xxx_types()` (i.e. `opflow.get_gen_bus_voltage_types()`). The `opflow.get_xxx_types` functions yield a list of the enum values.

Below are some code examples of how to get and use these values. 

**Code Examples**

Set with a string
```python
>>> opflow.set_initialization_type('OPFLOWINIT_FROMFILE')
>>> opflow.get_initialization_type()
OPFLOWInitializationType.OPFLOWINIT_FROMFILE
```

Set with an integer
```python
>>> opflow.set_initialization_type(1)
>>> opflow.get_initialization_type()
OPFLOWInitializationType.OPFLOWINIT_FROMFILE
```

Set with an enum instance
```python
>>> opflow.set_initialization_type(exago.OPFLOWInitializationType.OPFLOWINIT_FROMFILE)
>>> opflow.get_initialization_type()
OPFLOWInitializationType.OPFLOWINIT_FROMFILE
```

Set via getter function:
```python
>>> types = opflow.get_objective_types()
[<OPFLOWObjectiveType.MIN_GEN_COST: 0>, <OPFLOWObjectiveType.MIN_GENSETPOINT_DEVIATION: 1>, <OPFLOWObjectiveType.NO_OBJ: 2>]
>>> opflow.set_objective_type(types[0])
>>> opflow.get_objective_type()
OPFLOWObjectiveType.MIN_GEN_COST
```

See [the OPFLOW testing file](../tests/interfaces/python/test_opflow.py) for more usage examples. 

### Building

This documentation only applies to ExaGO >=v1.2.1.

When building ExaGO, enable the CMake option `EXAGO_ENABLE_PYTHON`.
For example,

```console
cmake .. \
  -D... # Set other ExaGO configuration options
  -DEXAGO_ENABLE_PYTHON=ON # Make sure Python bindings are enabled
make -j 12 install
```

This will install the Python bindings into the installation prefix in the
standard Python library suffix.

Using Python 3.8.0 on a Power9 system, ExaGO will install the Python bindings 
library `exago.cpython-38-powerpc64le-linux-gnu.so` into
`<install-prefix>/lib/python3.8/site-packages/exago.cpython-38-powerpc64le-linux-gnu.so`

## Environment

Currently ExaGO requires Python version greater than 3.6.
Once Python is available in your environment (`module load python`), you will need to add the ExaGO Python libraries to your `PYTHONPATH` according to your target Python version:

```bash
# For Python 3.6
export PYTHONPATH=/<exago_install>/lib/python3.6/site-packages:$PYTHONPATH
# For Python 3.7
export PYTHONPATH=/<exago_install>/lib/python3.7/site-packages:$PYTHONPATH
```

For ExaGO installations in non-standard locations, the `LD_LIBRARY_PATH` (`DYLD_LIBRARY_PATH` for MacOS) will also need to point to the ExaGO install location:

```bash
# For MacOS
export DYLD_LIBRARY_PATH=/<exago_install>/lib:$DYLD_LIBRARY_PATH
# For linux
export LD_LIBRARY_PATH=/<exago_install>/lib:$LD_LIBRARY_PATH
```

A sample opflow script [`test.py`](/interfaces/python/test.py) is provided, with `exago.prefix()` providing the path to the installation directory of ExaGO. For more example usage and for the tests that cover this code, see [`tests/interfaces/python`](/tests/interfaces/python).

## History

ExaGO pre v1.1 had optional Python bindings that could be enabled. HiOp, a
critical dependency, updated to Umpire v6 when GPU and RAJA options were
enabled. Because Umpire after v6 ships with CUDA device code, a final device
link step was required for any libraries/executables. This drastically
complicated the ctypes Python bindings, which relied on calling out to the
shared library directly.

ExaGO v1.2.1 was the first version to ship Python bindings that used Pybind11,
which creates a shared library directly importable from Python, which simplified
the user experience and made it possible to call ExaGO from Python when only 
static libraries are generated.
