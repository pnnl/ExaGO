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

### Dependency
To use this application, one must have ExaGO built with Ipopt. Even when using HiOp as the main solver, this application still uses Ipopt (to solve the base subproblem).

### Usage
The SOPFLOW code is executed via
```
mpiexec -n <N> ./sopflow <options>
```
where \<options\> are the available command line options as given in the next section.

### Options
The current version has several options available for SOPFLOW. These options can be set either through the options file `options/sopflowoptions` or via the command line.

|  Option Name | Description | Values (Default value) | Compatibility |
|:-----|:----|:-----|:-----|
|-netfile| Name of network file in MATPOWER format| ([case9mod.m](../../datafiles/case9/case9mod.m))|  4096 characters max. |
|-windgen| Name of wind scenario list file | ([10_scenarios_9bus.csv](../../datafiles/case9/10scenarios_9bus.csv)) | 4096 characters max. Uses a native format for describing scenarios. See [10_scenarios_9bus.csv](../../datafiles/case9/10scenarios_9bus.csv)|
|-ctgcfile| Name of contingency list file | |4096 characters max. Uses a native format for describing contingencies. |
|-sopflow_Ns | Number of scenarios || With this option set, SOPFLOW will only pick up the first Ns scenarios in the scenario file. To select all scenarios, use `Ns = -1` |
|-sopflow_model | SOPFLOW model type | GENRAMP, GENRAMPC (GENRAMP) | |
|-sopflow_solver | Optimization solver | (Ipopt), HiOp, or EMPAR | See the note below on solvers |
|-sopflow_mode | Mode of operation | 0 or 1 (0) | See the note below on mode of operation |
|-sopflow_subproblem_solver | Optimization solver for the subproblem when using HiOp solver| Ipopt or HiOp (Ipopt) | See [opflow](opflow.md) page for description of solvers |
|-sopflow_subproblem_model | Model for the subproblem when using HiOp solver| (POWER_BALANCE_POLAR) | See [opflow](opflow.md) page for available models |
|-sopflow_tolerance|Optimization solver tolerance | (1e-6) | All solvers |
|-scopflow_enable_multicontingency | Each scenario has multiple contingencies | 0 or 1 (0)| |
|-sopflow_Nc | Number of contingencies || With this option set, SOPFLOW will only pick up the first Nc contingencies in the contingency file. To select all contingencies, use `Nc = -1` |
|-sopflow_flatten_contingencies | Flattens out the scenario-contingency structure |0 or 1 (0)| Only used when multi-contingency is enabled |
|-print_output| Print SOPFLOW solution to screen| 0 or 1 (0)| All solvers |
|-save_output| Save SOPFLOW solution to file | 0 or 1 (0)| All solvers. Saves solution for each scenario. |

#### Scenarios
There are two types of uncertainties supported: Wind generation and load (currently not implemented).
Scenarios are specified in a native format. See `datafiles/case9/10scenarios_9bus.csv` as an example of a scenario file that describes wind generation scenarios for the 9-bus case.

#### Solver
SOPFLOW supports solving the problem using Ipopt, HiOP, or EMPAR solvers. With Ipopt, SOPFLOW can be only run on one processor (N = 1) as Ipopt only supports single process execution. HiOp supports solving the problem in parallel using a primal-decomposition algorithm. EMPAR is a parallel solver, however it merely executes an embarassingly parallel solver, i.e., all the scenarios are solved independently via optimal power flow.

#### Mode 
Set SOPFLOW to either run in `preventive` (0) or `corrective` (1) mode. In preventive mode, the base-case and scenario real-power dispatch is equal for the PV and PQ generators. Any power surplus/deficit is contributed by the swing generator only. In the corrective mode, the scenario real-power dispatch for all generators is allowed to deviate from the base-case limited by its 30-min ramping limit.
