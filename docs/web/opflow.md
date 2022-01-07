## AC optimal power flow (OPFLOW)
OPFLOW solves a nonlinear AC optimal power flow problem where the  (default) objective is generator-cost minimization subject to network balance, line flow, capacity, and voltage constaints. Load-loss and nodal power-imbalance can be (optionally) activated. 

### Usage

OPFLOW is executed via
```
mpiexec -n <N> ./opflow <options>
```
where \<options\> are the available command line options as given in the options section below. The opflow application is located in $EXAGO_INSTALL/bin where $EXAGO_INSTALL is the ExaGO install location. One can run OPFLOW from $EXAGO_INSTALL/bin or add $EXAGO_INSTALL/bin to the path for executing it from any other location.

*Note: If you are familiar with HPC tools, you may want to skip to the next section.*

In order to run OPFLOW in an HPC environment (such as PNNL cluster Newell), you
will have to request the required resources. For example, if OPFLOW has been compiled with support for GPUs, the user will have to request a GPU from the host job scheduler. On the PNNL cluster Newell, slurm is the host job scheduler, so an example command would be:

```console
$ srun -A <account name> -p newell --gres=gpu:1 -N 1 -n <N> ./opflow <options>
```

### Input

OPFLOW currently only supports input in the MATPOWER format. See [MATPOWER case format](https://matpower.org/docs/ref/matpower5.0/caseformat.html) for a description of the various fields for the network data. The `-netfile` option is used for specifying the input file. If no network file is specified as an input then OPFLOW defaults to running the 9-bus case with the network file `datafiles/case9/case9mod.m`.

```
./opflow -netfile <netfilename>
```

### Output

The output of OPFLOW can be either printed to screen (option `-print_output`) or saved to a file in MATPOWER format (option `-save_output`). With the `-save_output` option, the OPFLOW output is saved to a file named `opflowout.m` in the same directory from where the application is executed.

#### Print output
 ```
./opflow -netfile <netfilename> -print_output
```

#### Save output
 ```
./opflow -netfile <netfilename> -save_output
```

### Options

The behavior of OPFLOW is controlled through the different options given in the table below. If no options are provided, then the default options from [opflowoptions](../../options/opflowoptions) file are used. To see the available options, one can run OPFLOW with `-h` option and grep for opflow.
```
./opflow -h | grep opflow
```


|  Option Name | Description | Values (Default value) | Compatibility |
|:-----|:----|:-----|:-----|
|-opflow_model | Representation of network balance equations and bus voltages| See the different options for models below| |
|-opflow_solver | Optimization solver | See the different solver options below| |
|-opflow_initialization| Type of initialization| "MIDPOINT" (default)<br>"FROMFILE"<br>"ACPF"<br>"FLATSTART"| All solvers|
|-opflow_ignore_lineflow_constraints| Ignore line flow constraints| 0 or 1 (0)| All models|
|-opflow_include_loadloss_variables| Include load loss| 0 or 1 (0)| All models |
|-opflow_include_powerimbalance_variables| Allow power imbalance at buses| 0 or 1 (0)| All models |
|-opflow_loadloss_penalty| Penalty ($) for loss of load per load| (1000)| All models |
|-opflow_powerimbalance_penalty| Penalty ($) for  power imbalance at bus| (10000)| All models |
|-opflow_genbusvoltage| Control mode for generator bus voltage| "VARIABLE_WITHIN_BOUNDS"<br>"FIXED_WITHIN_QBOUNDS" (default)| POWER_BALANCE_POLAR model only |
|-opflow_has_gensetpoint| Real power set point set? | 0 or 1 (0)| All models |
|-opflow_objective| Objective function| "MIN_GEN_COST" (default)<br>"MIN_GENSETPOINT_DEVIATION"<br>"NO_OBJ"| MIN_GENSETPOINT_DEVIATION with POWER_BALANCE_POLAR model only |
|-opflow_use_agc| Use AGC for generator real power redispatch?| 0 or 1 (0)| POWER_BALANCE_POLAR model only |
|-opflow_tolerance|Optimization solver tolerance | (1e-6) | All solvers |
|-hiop_compute_mode|Controls the `-compute_mode` option for HIOP solver, i.e., where the HIOP solver computations run|"auto" (default)<br> "cpu"<br>"hybrid"<br>"gpu"| HiOp solver only |
|-hiop_verbosity_level|Controls the verbosity level for HIOP. 0 means no output is printed, 10 is max. output| 0 to 10 (0) See [HIOP verbosity levels](https://github.com/LLNL/hiop/blob/7e8adae9db757aed48e5c2bc448316307598258f/src/Utils/hiopLogger.hpp#L68)| HiOp solver only |
|-hiop_tolerance| HIOP solver tolerance| (1e-6)| HiOp solver only |
|-print_output| Print OPFLOW solution to screen| 0 or 1 (0)| All solvers |
|-save_output| Save OPFLOW solution to file | 0 or 1 (0)| All solvers |
|-netfile| Name of network file in MATPOWER format|4096 characters max. ([case9mod.m](../../datafiles/case9/case9mod.m))|

#### OPFLOW models

OPFLOW can be used with several different models i.e., representations of network balance equations and bus voltages. 

|  Model Name | Description |
|:-----|:----|
|POWER_BALANCE_POLAR|Power balance form with polar representation of voltages|
|POWER_BALANCE_HIOP| POWER_BALANCE_POLAR used only with HIOP solver on CPU|
|PBPOLRAJAHIOP| POWER_BALANCE_POLAR used only with HIOP solver on GPU|

#### OPFLOW optimization solvers

OPFLOW works with the following optimization solvers
1. [IPOPT](https://github.com/coin-or/Ipopt)
1. [HIOP](https://github.com/LLNL/hiop) - HiOp solver with mixed sparse-dense formulation
1. [HIOPSPARSE](https://github.com/LLNL/hiop) - HiOp solver with sparse formulation

IPOPT and HIOPSPARSE can run only on CPU, while HIOP can be run either on CPU or GPU.

#### Model-solver compatibility

When selecting a solver or a model, the following table should be referred to as not all solvers are compatible with all models.

|  Model name /Solver name | IPOPT | HIOP | HIOPSPARSE |
|:--------------------------|:-----:|:---------:|:-----:|
| POWER_BALANCE_POLAR  | x     |      |     x      |
| POWER_BALANCE_HIOP         |       | x  |   |
| PBPOLRAJAHIOP              |       | x    | |

#### Running OPFLOW on GPU
OPFLOW can execute on the GPU using the HIOP solver and PBPOLRAJAHIOP model. Note that both HIOP and ExaGO must be built with RAJA and UMPIRE libraries for executing on the GPU. There are two modes to execute OPFLOW on GPUs.

The `HYBRID` compute mode will run only the linear solver for the optimizer on GPU. All other calculations, including model evaluation, will be done on the CPU. 
```
mpiexec -n 1 ./opflow -opflow_solver HIOP -opflow_model PBPOLRAJAHIOP -hiop_compute_mode HYBRID <otheroptions>
```

The `GPU` compute mode will run _all_ calculations on the GPU.

```
mpiexec -n 1 ./opflow -opflow_solver HIOP -opflow_model PBPOLRAJAHIOP -hiop_compute_mode GPU <otheroptions>
```
