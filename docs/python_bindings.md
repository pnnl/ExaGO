## ExaGO Python Bindings

### Overview

The ExaGO Python bindings use an object-oriented API slightly different from the
C++ API.
The C++ API uses the application type in uppercase as the prefix for its
methods, where they are native methods in Python.

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

  /* Set up */
  ierr = OPFLOWSetUp(opflow);
  ExaGOCheckError(ierr);

  /* Solve */
  ierr = OPFLOWSolve(opflow);
  ExaGOCheckError(ierr);

  /* Save output to file 'opflowout' */
  ierr = OPFLOWSaveSolution(opflow, fmt, "opflowout");
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
exago.finalize()
```

Note that the Python bindings use snake case in accordance with standard Python
naming conventions.
Translating between the C++ API and the Python API should be relatively
straightforward.

***If you identify components of the C++ API that you need to call from Python,
please [open an issue on our issues page](https://gitlab.pnnl.gov/exasgd/frameworks/exago/-/issues).***

### Bindings Table

| C++ API | Python API | Notes |
|---|---|---|
| `ExaGOInitialize` | `exago.initialize` |  |
| `ExaGOFinalize` | `exago.finalize` |  |
| `PFLOW` | `exago.PFLOW` class | From here below, `pflow` objects are instances of the Python `exago.PFLOW` class. The same is the case for opflow, scopflow, tcopflow, and sopflow. |
| `PFLOWReadMatPowerData` | `pflow.read_mat_power_data` |  |
| `PFLOWSolve` | `pflow.solve` |  |
| `OPFLOW` | `exago.OPFLOW` class |  |
| `OPFLOWSetLoadLossPenalty` | `set_loadloss_penalty` | |
| `OPFLOWSetBusPowerImbalancePenalty` | `set_bus_powerimbalance_penalty` | |
| `OPFLOWSetTolerance` | `set_tolerance` | |
| `OPFLOWGetTolerance` | `get_tolerance` | |
| `OPFLOWGetObjective` | `get_objective` | |
| `OPFLOWSolve` | `solve` | |
| `OPFLOWPrintSolution` | `print_solution` | |
| `OPFLOWReadMatPowerData` | `read_mat_power_data` | |


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

