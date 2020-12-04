## Stochastic optimal power flow (SOPFLOW) -- IN DEVELOPMENT
SOPFLOW solves a stochastic optimal power flow problem. The problem is set up as a two-stage optimization problem where the first-stage (base-case) represents the normal operation of the grid (or the most likely forecast) and the second-stage comprises of $N_s$ scenarios of forecast deviation. Compactly, the problem can be set up in the following form:

```math
\begin{aligned}
\text{min}&~\sum_{i=0}^{N_s}p_if(x_i)& \\
&\text{s.t.}& \\
&~g(x_i) = 0~~~i \in \{0,N_s\}& \\
&~h(x_i) \le 0~~i \in \{0,N_s\}& \\
&x^- \le x_i \le x^+~~i\in \{0,N_s\}& \\
-\delta{x} & \le x_i - x_0 \le \delta{x}~~i \in \{1,N_s\}&
\end{aligned}
 ```

where $N_s$ is the number of scenarios and $p_i$ is the probability of the scenario. Each scenario is implemented as an optimal power flow formulation. See [OPFLOW](opflow.md). The last equation is the coupling between the 2nd stage forecast deviationscenarios and the base case. Depending on the `mode`, SOPFLOW can either be `preventive` (mode = 0) or `corrective` (mode = 1). In the preventive, the PV and PQ generator real power is not allowed to deviate from its base-case solution. The corrective mode allows deviation of the PV and PQ generator real power from the base-case dispatch constrained by its 30-min. ramp rate capability.


### Usage
The SOPFLOW code is executed via
```
mpiexec -n <N> ./sopflow <options>
```
where \<options\> are the available command line options as given in the next section.

### Options
The current version has several options available for SOPFLOW. These options can be either set through the options file `options/sopflowoptions` or via the command line.

#### Network file (-netfile <netfilename>): 
Set the name of the network file. Only MATPOWER format is currently supported.

```
mpiexec -n <N> ./sopflow -netfile <netfilename>
```

#### Scenario file (-scenfile <scenfilename>): 
Set the name of the file that describes the uncertainty scenarios. There are two types of uncertainties supported: Wind generation and load (currently not implemented). 
```
mpiexec -n <N> ./sopflow -netfile <netfilename> -scenfile <ctgcfilename>
```
Scenarios are specified in a native format. See `datafiles/TAM200_scenarios/scenarios.csv` as an example of a scenario file that describes wind generation scenarios for the TAMU 200-bus case.

#### Solver (-sopflow_solver <IPOPT or EMPAR>)
Set the solver to be used for SOPFLOW. With IPOPT, SOPFLOW can be only run on one processor (N = 1) as IPOPT only supports single process execution. EMPAR is a parallel solver, however it merely executes an embarassingly parallel solver, i.e., all the scenarios are solved independently via optimal power flow.
```
mpiexec -n <N> ./sopflow -netfile <netfilename> -ctgcfile <ctgcfilename> -sopflow_solver <IPOPT>
```
#### Mode (-sopflow_mode <0 or 1>)
Set SOPFLOW to either run in `preventive` (0) or `corrective` (1) mode. In preventive mode, the base-case and scenario real-power dispatch is equal for the PV and PQ generators. In the corrective mode, the scenario real-power dispatch for these generators is allowed to deviate from the base-case limited by its 30-min ramping limit. 
```
mpiexec -n <N> ./sopflow -netfile <netfilename> -scenfile <scenfilename> -sopflow_mode <0 or 1>
```

#### Number of contingencies (-sopflow_Ns <Ns>): 
Sets the number of scenarios. This should be less than or equal to the number of scenarios set in the scenario file.

```
./sopflow -netfile <netfilename> -scenfile <scenfilename> -sopflow_Ns <Ns>
```

With this option set, SOPFLOW will only pick up the first Ns scenarios in the scenario file. To select all scenarios, use `Ns = -1`

#### OPFLOW model (-opflow_model <modelname>)
This is an option inherited from [OPFLOW](opflow.md). It sets the model (representation of variables and equations) to be used for SOPFLOW. The default formulation is power balance form with polar representation for voltages
```
mpiexec -n <N> ./sopflow -netfile <netfilename> -scenfile <scenfilename> -opflow_model POWER_BALANCE_POLAR
```

All other OPFLOW options can be used with SOPFLOW to select the appropriate model and parameter options.
