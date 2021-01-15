## AC optimal power flow (OPFLOW)
OPFLOW solves a nonlinear AC optimal power flow problem. The  objective of this OPFLOW problem is generator-cost minimization subject to network balance, line flow, capacity, and voltage constaints. Load-loss and nodal power-imbalance can be (optionally) activated. OPFLOW can be used with several different models i.e., representations of network balance equations and bus voltages.
1. [**POWER_BALANCE_POLAR**](../opflow/pbpol.md): Power balance form with polar representation of voltages.
1. [**POWER_BALANCE_CARTESIAN**](../opflow/pbcar.md): Power balance form with cartesian representation of voltages.
1. [**CURRENT_BALANCE_CARTESIAN**]((../opflow/ibcar.md)): Current balance form with cartesian representation of voltages.
1. **CURRENT_BALANCE_CARTESIAN2**: Modified current balance form with bounds on line currents instead of MVA flow.
1. **POWER_BALANCE_HIOP**: An alternate version of power balance polar form used only with HIOP solver.
1. **PBPOLRAJAHIOP**: An alternate version of power balance polar for GPU rewritten with RAJA library used only with HIOP solver.

OPFLOW works with the following optimization solvers.
1. [IPOPT](https://github.com/coin-or/Ipopt)
1. [HiOp](https://github.com/LLNL/hiop) 
1. [PETSc/TAO](https://www.mcs.anl.gov/petsc/)

The compatibility of the different models and solvers is given in the table below:


|  Model name /Solver name             | IPOPT | HIOP | TAO | 
|:--------------------------:|:-----:|:----:|:---------:|
| POWER_BALANCE_POLAR        | x     |      |           |
| POWER_BALANCE_CARTESIAN    | x     |      | x         |
| CURRENT_BALANCE_CARTESIAN  | x     |      |           |
| CURRENT_BALANCE_CARTESIAN2 | x     |      |           |
| POWER_BALANCE_HIOP         |       | x    |           |
| PBPOLRAJAHIOP              |       | x    | x         |



### Usage
OPFLOW is executed via
```
mpiexec -n <N> ./opflow <options>
```
where \<options\> are the available command line options as given in the following sections. The opflow application is located in $EXAGO_INSTALL/bin where $EXAGO_INSTALL is the ExaGO install location.

Note that only TAO supports parallel execution. For all other solvers, use N = 1. 

#### Using GPUs
OPFLOW can execute on the GPU using the HIOP solver and PBPOLRAJAHIOP model. Note that both HIOP and ExaGO must be built with RAJA and UMPIRE libraries for executing on the GPU. The execution command for running OPFLOW on GPU is
```
mpiexec -n 1 ./opflow -opflow_solver HIOP -opflow_model PBPOLRAJAHIOP -hiop_compute_mode HYBRID <otheroptions>
```

### Run-time options
OPFLOW has several run-time options. These options can be either set through the options file `options/opflowoptions` or via the command line. To see the different options, run opflow with the help flag.
```
mpiexec -n <N> ./opflow -h | grep opflow
```
This will list the different options available.

#### Network file (-netfile <netfilename>): 
Set the name of the network file (MATPOWER format only currently). There is support for reading the PSSE format as well, but the PSSE raw data file does not contain the generator cost data. A separate file needs to be set for the generator cost which is not supported yet with OPFLOW.

```
mpiexec -n <N> ./opflow -netfile <netfilename>
```

#### Solver (-opflow_solver <IPOPT, TAO, HIOP>)
Set the solver to be used for OPFLOW. OPFLOW currently supports three solvers - IPOPT, TAO, and HiOp. If IPOPT or HiOp is chosen then OPFLOW can be only run on one processor (N = 1) as both these solvers only support single process execution. TAO supports parallel execution. The above table should be refered when choosing the appropriate solver and model.
```
mpiexec -n <N> ./opflow -netfile <netfilename> -opflow_solver <IPOPT, TAO, HiOp>
```
The default solver is IPOPT. If IPOPT is not used then the solver defaults to TAO

#### Set the model (-opflow_model <modelname>): 
Sets the model (represent of network balance equations and voltages). See the above table for the different models available.
```
./opflow -netfile <netfilename> -opflow_model <modelname>
```

#### Set the initialization type (-opflow_initialization <initializationtype>)
Sets the initialization (initial guess) for OPFLOW.
```
./opflow -netfile <netfilename> -opflow_initialization <initializationtype>
```
The initializations supported are
- `MIDPOINT` - Uses mid-point of the bounds for voltages and generator powers
- `FLATSTART` - Flat start (1.0) for voltages, mid-point for generator powers
- `ACPF` - Initialization using AC power flow
- `FROMFILE` - Uses values given in file

#### Ignore line flow constraints (-opflow_ignore_lineflow_constraints <0,1>)
With this option set, OPFLOW ignores all line flow constraints
```
./opflow -netfile <netfilename> -opflow_ignore_lineflow_constraints
```

#### Include load loss (-opflow_include_loadloss_variables <0,1>)
With this option set, OPFLOW can also calculate the loss of load. Extra variables for calculating the load loss at each bus are inserted in OPFLOW and the objective function is modified to include the load loss.
```
./opflow -netfile <netfilename> -opflow_include_loadloss_variables
```

#### Load loss penalty parameter(-opflow_loadloss_penalty <penalty_cost>)
Through this option, the penalty cost for load loss when `-opflow_include_loadloss_variables` can be set. The default is 1000. 
```
./opflow -netfile <netfilename> -opflow_include_loadloss_variables -opflow_loadlloss_penalty <penalty_cost>
```

#### Include power imbalance variables (-opflow_include_powerimbalance_variables <0,1>)
With this option set, OPFLOW inserts extra variables for measuring the power imbalance (real and reactive power) at buses. For a feasible optimal power flow solution, the power imbalance at each bus is 0. In case there is infeasibity, the power imbalance variables measure the amount of infeasability at each bus. 
```
./opflow -netfile <netfilename> -opflow_include_powerimbalance_variables
```

#### Power imbalance penalty parameter(-opflow_powerimbalance_penalty <penalty_cost>)
Through this option, the penalty cost for power imbalance when `-opflow_include_powerimbalance_variables` can be set. The default is 1000. 
```
./opflow -netfile <netfilename> -opflow_include_powerimbalance_variables -opflow_powerimbalance_penalty <penalty_cost>
```

#### Tolerance (-opflow_tolerance <tolerance>)
This option sets the optimization tolerance. Default is 1e-6.