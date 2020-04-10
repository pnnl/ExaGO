## Security-constrained optimal power flow (SCOPFLOW)
SCOPFLOW solves a contingency-constrained optimal power flow problem. The problem is set up as a two-stage optimization problem where the first-stage (base-case) represents the normal operation of the grid and the second-stage comprises of $N_s$ scenarios representing deviations from the normal operation such as those caused by contingencies or renewable uncertainties. Compactly, the problem can be set up in the following form:

```math
\begin{aligned}
\text{min}&~\sum_{i=0}^{Ns} f(x_i)& \\
&\text{s.t.}& \\
&~g(x_i) = 0~~~i \in \{0,N_s\}& \\
&~h(x_i) \le 0~~i \in \{0,N_s\}& \\
&x^- \le x_i \le x^+~~i\in \{0,N_s\}& \\
-\delta{x}& \le x_i - x_0 \le \delta{x}~~i \in \{1,N_s\}&
\end{aligned}
 ```
where $N_s$ is the number of scenarios. Each scenario is an optimal power flow formulation. See [OPFLOW](opflow.md). The last equation is the coupling between the 2nd stage scenarios and the first-stage that enforces a limit on the deviation on second stage decision variable $x_i$ from its corresponding base case decision variable $x_0$.

### Usage
The SCOPFLOW code is executed via
```
mpiexec -n <N> ./SCOPFLOW <options>
```
where \<options\> are the available command line options as given in the next section. 

### Options
The current version has several options available for SCOPFLOW. These options can be either set through the options file `options/scopflowoptions` or via the command line.

#### Network file (-netfile <netfilename>): 
Set the name of the network file (MATPOWER format only currently). There is support for reading the PSSE format as well, but the PSSE raw data file does not contain the generator cost data. A separate file needs to be set for the generator cost which is not supported yet with SCOPFLOW.

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

#### Number of scenarios (-scopflow_Ns <Ns>): 
Sets the number of second-stage scenarios. This should be less than or equal to the number of contingencies set in the contingency file.
```
./SCOPFLOW -netfile <netfilename> -ctgcfile <ctgcfilename> -scopflow_Ns <Ns>
```
With this option set, SCOPFLOW will only pick up the first Ns contingencies in the contingency file. 
With PIPS as the solver for SCOPFLOW, the number of scenarios set should be larger than the number of ranks ($N_s > N$).

#### Include second stage cost (-scopflow_first_stage_gen_cost_only <0,1>)
If this option is set then only the first stage (base case) objective is considered, second stage (contingency scenarios) objective costs are ignored. 
```
mpiexec -n <N> ./SCOPFLOW -netfile <netfilename> -ctgcfile <ctgcfilename> -scopflow_first_stage_gen_cost_only <0,1>
```

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





