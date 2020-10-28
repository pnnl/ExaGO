## Security-constrained optimal power flow (SCOPFLOW)
SCOPFLOW solves a contingency-constrained optimal power flow problem. The problem is set up as a two-stage optimization problem where the first-stage (base-case) represents the normal operation of the grid and the second-stage comprises of $N_c$ contingency scenarios. Compactly, the problem can be set up in the following form:

```math
\begin{aligned}
\text{min}&~f(x_0)& \\
&\text{s.t.}& \\
&~g(x_i) = 0~~~i \in \{0,N_c\}& \\
&~h(x_i) \le 0~~i \in \{0,N_c\}& \\
&x^- \le x_i \le x^+~~i\in \{0,N_c\}& \\
-\delta{x} & \le x_i - x_0 \le \delta{x}~~i \in \{1,N_c\}&
\end{aligned}
 ```
where $N_c$ is the number of contingency. Each scenario is an optimal power flow formulation. See [OPFLOW](opflow.md). The last equation is the coupling between the 2nd stage contingency scenarios and the first-stage. Depending on the `mode`, SCOPFLOW can either be `preventive` (mode = 0) or `corrective` (mode = 1). In the preventive, the PV and PQ generator real power is not allowed to deviate from its base-case solution. The corrective mode allows deviation of the PV and PQ generator real power from the base-case dispatch constrained by its 30-min. ramp rate capability.


### Usage
The SCOPFLOW code is executed via
```
mpiexec -n <N> ./scopflow <options>
```
where \<options\> are the available command line options as given in the next section.

### Options
The current version has several options available for SCOPFLOW. These options can be either set through the options file `options/scopflowoptions` or via the command line.

#### Network file (-netfile <netfilename>): 
Set the name of the network file. Only MATPOWER format is currently supported.

```
mpiexec -n <N> ./scopflow -netfile <netfilename>
```

#### Contingency file (-ctgcfile <ctgcfilename>): 
Set the name of the contingency data file. The contingency scenarios are set through this file.
```
mpiexec -n <N> ./scopflow -netfile <netfilename> -ctgcfile <ctgcfilename>
```
Contingencies can either be specified in PTI format (.con file) or a native format. The description of the native format is given in the header file `include/scopflow.h`. SCOPFLOW supports single/multiple generator and line/transformer outage contingencies.

#### Solver (-scopflow_solver <IPOPT>)
Set the solver to be used for SCOPFLOW. Currently, only IPOPT is supported. With IPOPT, SCOPFLOW can be only run on one processor (N = 1) as IPOPT only supports single process execution.
```
mpiexec -n <N> ./scopflow -netfile <netfilename> -ctgcfile <ctgcfilename> -scopflow_solver <IPOPT>
```

#### Mode (-scopflow_mode <0 or 1>)
Set SCOPFLOW to either run in `preventive` (0) or `corrective` (1) mode. In preventive mode, the base-case and contingency real-power dispatch is equal for the PV and PQ generators. In the corrective mode, the contingency real-power dispatch for these generators is allowed to deviate from the base-case limited by its 30-min ramping limit. 
```
mpiexec -n <N> ./scopflow -netfile <netfilename> -ctgcfile <ctgcfilename> -scopflow_mode <0 or 1>
```

#### Number of contingencies (-scopflow_Nc <Nc>): 
Sets the number of contingencies. This should be less than or equal to the number of contingencies set in the contingency file.
```
./scopflow -netfile <netfilename> -ctgcfile <ctgcfilename> -scopflow_Nc <Nc>
```
With this option set, SCOPFLOW will only pick up the first Nc contingencies in the contingency file. To select all contingencies, use `Nc = -1`

#### OPFLOW model (-opflow_model <modelname>)
This is an option inherited from [OPFLOW](opflow.md). It sets the model (representation of variables and equations) to be used for SCOPFLOW. The default formulation is power balance form with polar representation for voltages
```
mpiexec -n <N> ./scopflow -netfile <netfilename> -ctgcfile <ctgcfilename> -opflow_model POWER_BALANCE_POLAR
```

All other OPFLOW options can be used with SCOPFLOW to select the appropriate model and parameter options.
