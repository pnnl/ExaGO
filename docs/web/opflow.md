## AC optimal power flow (OPFLOW)
OPFLOW solves a nonlinear AC optimal power flow problem. The  objective of this OPFLOW problem is generator-cost minimization subject to netwok balance,line flow,  capacity,and voltage constaints. Load-loss and nodal power-imbalance can be (optionally) activated. OPFLOW can be used with four different representations of network balance equations and bus voltages.
1. [POWER_BALANCE_POLAR](../opflow/pbpol.md): Power balance form with polar representation of voltages.
1. [POWER_BALANCE_CARTESIAN](../opflow/pbcar.md): Power balance form with cartesian representation of voltages.
1. [CURRENT_BALANCE_CARTESIAN]((../opflow/ibcar.md)): Current balance form with cartesian representation of voltages.
1. CURRENT_BALANCE_CARTESIAN2: Modified current balance form with bounds on line currents instead of MVA flow.

### Usage
OPFLOW is executed via
```
mpiexec -n <N> ./OPFLOW <options>
```
where \<options\> are the available command line options as given in the next section.

Note that only TAO supports parallel execution. For all other solvers, use N = 1. 

### Options
The current version has several options available for OPFLOW. These options can be either set through the options file `options/opflowoptions` or via the command line. To see the different options, run OPFLOW with the help flag.
```
mpiexec -n <N> ./OPFLOW -h | grep OPFLOW
```
This will list the different options available.

#### Network file (-netfile <netfilename>): 
Set the name of the network file (MATPOWER format only currently). There is support for reading the PSSE format as well, but the PSSE raw data file does not contain the generator cost data. A separate file needs to be set for the generator cost which is not supported yet with OPFLOW.

```
mpiexec -n <N> ./OPFLOW -netfile <netfilename>
```

#### Solver (-opflow_solver <IPOPT, TAO, HIOP>)
Set the solver to be used for OPFLOW. OPFLOW currently supports three solvers - IPOPT, TAO, and HiOp. If IPOPT or HiOp is chosen then OPFLOW can be only run on one processor (N = 1) as both these solvers only support single process execution. TAO supports parallel execution.
```
mpiexec -n <N> ./OPFLOW -netfile <netfilename> -scopflow_solver <IPOPT, TAO, HiOp>
```
The default solver is TAO.

#### Set the formulation (-opflow_formulation <formulationname>): 
Sets the formulation (represent of network balance equations and voltages). There are four different representations that can be run as mentioned above
```
./OPFLOW -netfile <netfilename> -opflow_formulation <formulationname>
```
The formulation name can be one of `POWER_BALANCE_POLAR`,`POWER_BALANCE_CARTESIAN`,`CURRENT_BALANCE_CARTESIAN`, or `CURRENT_BALANCE_CARTESIAN2`. The default is `POWER_BALANCE_CARTESIAN`.

#### Set the initialization type (-opflow_initialization <initializationtype>)
Sets the initialization (initial guess) for OPFLOW.
```
./OPFLOW -netfile <netfilename> -opflow_initialization <initializationtype>
```
The initializations supported are
- `MIDPOINT` - Uses mid-point of the bounds for voltages and generator powers
- `FLATSTART` - Flat start (1.0) for voltages, mid-point for generator powers
- `ACPF` - Initialization using AC power flow
- `FROMFILE` - Uses values given in file

#### Ignore line flow constraints (-opflow_ignore_lineflow_constraints <0,1>)
With this option set, OPFLOW ignores all line flow constraints
```
./OPFLOW -netfile <netfilename> -opflow_ignore_lineflow_constraints
```

#### Include load loss (-opflow_include_loadloss_variables <0,1>)
With this option set, OPFLOW can also calculate the loss of load. Extra variables for calculating the load loss at each bus are inserted in OPFLOW and the objective function is modified to include the load loss.
```
./OPFLOW -netfile <netfilename> -opflow_include_loadloss_variables
```

#### Load loss penalty parameter(-opflow_loadloss_penalty <penalty_cost>)
Through this option, the penalty cost for load loss when `-opflow_include_loadloss_variables` can be set. The default is 1000. 
```
./OPFLOW -netfile <netfilename> -opflow_include_loadloss_variables -opflow_loadlloss_penalty <penalty_cost>
```

#### Include power imbalance variables (-opflow_include_powerimbalance_variables <0,1>)
With this option set, OPFLOW inserts extra variables for measuring the power imbalance (real and reactive power) at buses. For a feasible optimal power flow solution, the power imbalance at each bus is 0. In case there is infeasibity, the power imbalance variables measure the amount of infeasability at each bus. 
```
./OPFLOW -netfile <netfilename> -opflow_include_powerimbalance_variables
```

#### Power imbalance penalty parameter(-opflow_powerimbalance_penalty <penalty_cost>)
Through this option, the penalty cost for power imbalance when `-opflow_include_powerimbalance_variables` can be set. The default is 1000. 
```
./OPFLOW -netfile <netfilename> -opflow_include_powerimbalance_variables -opflow_powerimbalance_penalty <penalty_cost>
```
