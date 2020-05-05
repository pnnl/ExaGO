## Security-constrained optimal power flow (SCOPFLOW)
SCOPFLOW solves a contingency-constrained optimal power flow problem. The problem is set up as a two-stage optimization problem where the first-stage (base-case) represents the normal operation of the grid and the second-stage comprises of $N_c$ contingency scenarios. Compactly, the problem can be set up in the following form:

```math
\begin{aligned}
\text{min}&~f(x_0)& \\
&\text{s.t.}& \\
&~g(x_i) = 0~~~i \in \{0,N_c\}& \\
&~h(x_i) \le 0~~i \in \{0,N_c\}& \\
&x^- \le x_i \le x^+~~i\in \{0,N_c\}& \\
-\text{mode}*\delta{x} & \le x_i - x_0 \le \text{mode}*\delta{x}~~i \in \{1,N_c\}&
\end{aligned}
 ```
where $N_c$ is the number of contingency. Each scenario is an optimal power flow formulation. See [OPFLOW](opflow.md). The last equation is the coupling between the 2nd stage contingency scenarios and the first-stage. Depending on the `mode`, SCOPFLOW can either be `preventive` (mode = 0) or `corrective` (mode = 1). In the preventive, the generator real power is not allowed to deviate from its base-case solution. The corrective mode allows deviation from the base-case dispatch constrained its 30-min. ramp rate capability.


### Usage
The SCOPFLOW code is executed via
```
mpiexec -n <N> ./SCOPFLOW <options>
```
where \<options\> are the available command line options as given in the next section.

### Options
The current version has several options available for SCOPFLOW. These options can be either set through the options file `options/scopflowoptions` or via the command line.

#### Network file (-netfile <netfilename>): 
Set the name of the network file. Only MATPOWER format is currently supported.

```
mpiexec -n <N> ./SCOPFLOW -netfile <netfilename>
```

#### Contingency file (-ctgcfile <ctgcfilename>): 
Set the name of the contingency data file. The contingency scenarios are set through this file.
```
mpiexec -n <N> ./SCOPFLOW -netfile <netfilename> -ctgcfile <ctgcfilename>
```
The format of the contingencies is explained in the header file `include/scopflow.h`

#### Solver (-scopflow_solver <IPOPT or PIPS>)
Set the solver to be used for SCOPFLOW. Currently, only IPOPT and PIPS-NLP are supported. If IPOPT is chosen then SCOPFLOW can be only run on one processor (N = 1) as IPOPT only supports single process execution. PIPS supports parallel execution.
```
mpiexec -n <N> ./SCOPFLOW -netfile <netfilename> -ctgcfile <ctgcfilename> -scopflow_solver <IPOPT or PIPS>
```
If this option is not set then SCOPFLOW uses IPOPT solver.

Note: Currently, only IPOPT solver is supported.

#### Mode (-scopflow_mode <0 or 1>)
Set SCOPFLOW to either run in `preventive` (0) or `corrective` (1) mode. In preventive mode, the base-case and contingency real-power dispatch is equal. In the corrective mode, the contingency real-power dispatch is allowed to deviate from the base-case maximum upto its 30-min ramping limit. 
```
mpiexec -n <N> ./SCOPFLOW -netfile <netfilename> -ctgcfile <ctgcfilename> -scopflow_mode <0 or 1>
```

#### Number of scenarios (-scopflow_Nc <Nc>): 
Sets the number of second-stage scenarios. This should be less than or equal to the number of contingencies set in the contingency file.
```
./SCOPFLOW -netfile <netfilename> -ctgcfile <ctgcfilename> -scopflow_Nc <Nc>
```
With this option set, SCOPFLOW will only pick up the first Nc contingencies in the contingency file. To select all contingencies, use `Nc = -1` 
With PIPS as the solver for SCOPFLOW, the number of scenarios set should be larger than the number of ranks ($N_c > N$).

#### Formulation (-opflow_formulation <formulationname>)
This is an option inherited from [OPFLOW](opflow.md). It sets the formulation (representation of variables and equations) to be used for SCOPFLOW. The default formulation is power balance form with polar representation for voltages
```
mpiexec -n <N> ./SCOPFLOW -netfile <netfilename> -ctgcfile <ctgcfilename> -opflow_formulation <formulationname>
```
Currently, the following formulations are supported by SCOPFLOW.
1. [POWER_BALANCE_POLAR](../opflow/pbpol.md): Power balance form with polar representation of voltages.
1. [POWER_BALANCE_CARTESIAN](../opflow/pbcar.md): Power balance form with cartesian representation of voltages.
1. [CURRENT_BALANCE_CARTESIAN]((../opflow/ibcar.md)): Current balance form with cartesian representation of voltages.
1. CURRENT_BALANCE_CARTESIAN2: Modified current balance form with bounds on line currents instead of MVA flow.





