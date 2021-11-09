## Stochastic security-constrained multi-period optimal power flow (SOPFLOW)
SOPFLOW solves a stochastic security-constrained multi-period optimal power flow problem. The problem is set up as a two-stage optimization problem where the first-stage (base-case) represents the normal operation of the grid (or the most likely forecast) and the second-stage comprises of $`N_s`$ scenarios of forecast deviation. Compactly, the problem can be set up in the following form:

```math
\begin{aligned}
\text{min}&~\sum_{s=0}^{N_s-1}\pi_s\sum_{c=0}^{N_c-1}\sum_{t=0}^{N_t-1}f(x_{s,c,t})& \\
&\text{s.t.}& \\
&~g(x_{s,c,t}) = 0,                                        &s \in \{1,N_s-1\},c \in \{0,N_c-1\}, t \in \{0,N_t-1\}& \\
&~h(x_{s,c,t}) \le 0,                                      &s \in \{1,N_s-1\},c \in \{0,N_c-1\}, t \in \{0,N_t-1\}& \\
x^- & \le x_{s,c,t} \le x^+,                               &s \in \{1,N_s-1\},c \in \{0,N_c-1\}, t \in \{0,N_t-1\}& \\
-\Delta x_t & \le x_{s,c,t} - x_{s,c,t-\Delta{t}} \le \Delta x_t,&s \in \{1,N_s-1\},c \in \{0,N_c-1\}, t \in \{1,N_t-1\}&\\
-\Delta x_c & \le x_{s,c,0} - x_{s,0,0} \le \Delta x_c,&s \in \{1,N_s-1\},c \in \{1,N_c-1\}& \\
-\Delta x_s & \le x_{s,0,0} - x_{0,0,0} \le \Delta x_s,&s \in \{1,N_s-1\}&
\end{aligned}
```

where $`N_s`$ is the number of scenarios, $`N_c`$ is the number of contingencies, $`N_t`$ is the number of periods and $`\pi_s`$ is the probability of the scenario. Each scenario can either be an optimal power flow (see [OPFLOW](opflow.md)) or security-constrained optimal power flow (see [SCOPFLOW](scopflow.md)). Each security-constrained optimal power flow case can be single or multi-period. The last equation is the coupling between the 2nd stage forecast deviations and the base case. Depending on the `mode`, SOPFLOW can either be `preventive` (mode = 0) or `corrective` (mode = 1). In the preventive mode, generator real power output is fixed to the base-case values except for any renewable generation (wind, solar), and generators at reference bus(es).The corrective mode allows deviation of the PV and PQ generator real power from the base-case dispatch constrained by its 30-min. ramp rate capability.

Depending on the options chosen, SOPFLOW can be run in three modes:
1. Single-period no-contingency stochastic optimal power flow
1. Single-period multi-contingency stochastic optimal power flow
1. Multi-period multi-contingency stochastic optimal power flow

### Usage
The SOPFLOW code is executed via
```
mpiexec -n <N> ./sopflow <options>
```
where \<options\> are the available command line options as given in the next section.

### Options
The current version has several options available for SOPFLOW. These options can be either set through the options file `options/sopflowoptions` or via the command line.

#### Network file (-netfile \<netfilename\>): 
Set the name of the network file. Only MATPOWER format is currently supported. 4096 characters max.

```
mpiexec -n <N> ./sopflow -netfile <netfilename>
```

#### Scenario file (-scenfile \<scenfilename\>): 
Set the name of the file that describes the uncertainty scenarios. There are two types of uncertainties supported: Wind generation and load (currently not implemented). 4096 characters max.
```
mpiexec -n <N> ./sopflow -netfile <netfilename> -scenfile <ctgcfilename>
```
Scenarios are specified in a native format. See `datafiles/TAM200_scenarios/scenarios.csv` as an example of a scenario file that describes wind generation scenarios for the TAMU 200-bus case.

#### Solver (-sopflow_solver \<IPOPT or EMPAR\>)
Set the solver to be used for SOPFLOW. With IPOPT, SOPFLOW can be only run on one processor (N = 1) as IPOPT only supports single process execution. EMPAR is a parallel solver, however it merely executes an embarassingly parallel solver, i.e., all the scenarios are solved independently via optimal power flow.
```
mpiexec -n <N> ./sopflow -netfile <netfilename> -ctgcfile <ctgcfilename> -sopflow_solver <IPOPT>
```
#### Mode (-sopflow_mode \<0 or 1\>)
Set SOPFLOW to either run in `preventive` (0) or `corrective` (1) mode. In preventive mode, the base-case and scenario real-power dispatch is equal for the PV and PQ generators. In the corrective mode, the scenario real-power dispatch for these generators is allowed to deviate from the base-case limited by its 30-min ramping limit. 
```
mpiexec -n <N> ./sopflow -netfile <netfilename> -scenfile <scenfilename> -sopflow_mode <0 or 1>
```

#### Number of scenarios (-sopflow_Ns \<Ns\>): 
Sets the number of scenarios. This should be less than or equal to the number of scenarios set in the scenario file.

```
./sopflow -netfile <netfilename> -scenfile <scenfilename> -sopflow_Ns <Ns>
```

With this option set, SOPFLOW will only pick up the first Ns scenarios in the scenario file. To select all scenarios, use `Ns = -1`

#### Enable multi-contingency [Only for multi-contingency SOPFLOW] (-sopflow_enable_multicontingency \<0 or 1\>)
Disable/Enable multicontingency SOPFLOW (disabled by default). All SCOPFLOW options (see [SCOPFLOW](scopflow.md)) can be used when the multi-contingency option is enabled.

#### Enable multi-period SCOPFLOW (-scopflow_enable_multiperiod \<0 or 1\>)
Disable/Enable multi-period multicontingency SCOPFLOW (disabled by default). All options for multi-period SCOPFLOW (see [SCOPFLOW](scopflow.md)) can be used when the multi-period option is enabled.
