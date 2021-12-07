## Security-constrained optimal power flow (SCOPFLOW)
SCOPFLOW solves a contingency-constrained optimal power flow problem. The problem is set up as a two-stage optimization problem where the first-stage (base-case) represents the normal operation of the grid and the second-stage comprises of $`N_c`$ contingency scenarios. Compactly, the problem can be set up in the following form:

```math
\begin{aligned}
\text{min}&~\sum_{c \in \{0,N_c\}}f_c(x_c)& \\
&\text{s.t.}& \\
&~g_c(x_c) = 0~~~c \in \{0,N_c\}& \\
&~h_c(x_c) \le 0~~c \in \{0,N_c\}& \\
&x^- \le x_c \le x^+~~c\in \{0,N_c\}& \\
-\Delta{x}_c & \le x_c - x_0 \le \Delta{x}_c~~c \in \{1,N_c\}&
\end{aligned}
```

where $`N_c`$ is the number of contingencies. The total number of scenarios equals $`N_c + 1`$, i.e., the base-case + $`N_c`$ contingencies. Each scenario is an optimal power flow formulation. See [OPFLOW](opflow.md). The last equation is the coupling between the 2nd stage contingency scenarios and the first-stage.  Each contingency scenario can either be single-period or multi-period. In the multi-period mode, additional data files for the load and wind generation profiles can be set via command line options. Multi-period SCOPFLOW is activated either by setting the command line option `-scopflow_enable_multiperiod` OR calling `SCOPFLOWEnableMultiPeriod`.

Depending on the `mode`, SCOPFLOW can either be `preventive` (mode = 0) or `corrective` (mode = 1). In the preventive mode, the PV and PQ generator real power is fixed to its correspoinding base-case values. Any power offset/make-up is done by the swing bus generators. The corrective mode allows deviation of the PV and PQ generator real power from the base-case dispatch, constrained by its 30-min. ramp rate capability.


### Usage
SCOPFLOW is executed via
```
mpiexec -n <N> ./scopflow <options>
```
where \<options\> are the available command line options as given in the next section.

### Options
The current version has several options available for SCOPFLOW. These options can be set either through the options file `options/scopflowoptions` or via the command line.

|  Option Name | Description | Values (Default value) | Compatibility |
|:-----|:----|:-----|:-----|
|-netfile| Name of network file in MATPOWER format| ([case9mod.m](../../datafiles/case9/case9mod.m))|  4096 characters max. |
|-ctgcfile| Name of contingency list file | ([case9.cont](../../datafiles/case9/case9.cont)) | 4096 characters max. Uses a native format for describing contingencies. See [scopflow.h](../../include/scopflow.h)|
|-scopflow_Nc | Number of contingencies || With this option set, SCOPFLOW will only pick up the first Nc contingencies in the contingency file. To select all contingencies, use `Nc = -1` |
|-scopflow_solver | Optimization solver | (IPOPT), HIOP, or EMPAR | See the note below on solvers |
|-scopflow_mode | Mode of operation | 0 or 1 (0) | See the note below on mode of operation |
|-scopflow_subproblem_solver | Optimization solver for the subproblem when using HIOP solver| IPOPT or HIOP (IPOPT) | See [opflow](opflow.md) page for description of solvers |
|-scopflow_subproblem_model | Model for the subproblem when using HIOP solver| (POWER_BALANCE_POLAR) | See [opflow](opflow.md) page for available models |
|-scopflow_mode | Mode of operation | 0 or 1 (0) | See the note below on mode of operation |
|-scopflow_tolerance|Optimization solver tolerance | (1e-6) | All solvers |
|-scopflow_enable_multiperiod | Include multi-period | 0 or 1 (0)| Only compatible with IPOPT solver |
|-scopflow_duration| Duration for multi-period run in hours| | Only when multi-period is enabled |
|-scopflow_dt| Time-step for multi-period run in minutes| | Only when multi-period is enabled |
|-scopflow_ploadprofile| Active power load profile filename| | Only when multi-period is enabled |
|-scopflow_qloadprofile| Reactive power load profile filename| | Only when multi-period is enabled |
|-scopflow_windgenprofile| Wind generation profile filename| | Only when multi-period is enabled |
|-print_output| Print SCOPFLOW solution to screen| 0 or 1 (0)| All solvers |
|-save_output| Save OPFLOW solution to file | 0 or 1 (0)| All solvers |

#### Contingencies 
Contingencies can either be specified in PTI format (.con file) or a native format. The description of the native format is given in the header file `include/scopflow.h`. SCOPFLOW supports single/multiple generator and line/transformer outage contingencies.

#### Solver
SCOPFLOW can be solved with either IPOPT, HIOP, or EMPAR. With IPOPT, SCOPFLOW can be only run on one processor (N = 1) as IPOPT only supports single process execution. HIOP supports a distributed solution allowing SCOPFLOW to be solved in parallel. It uses a two-stage primal decomposition algorithm for solving the problem. 

In addition, one can solve SCOPFLOW in an embarrasingly parallel model with the EMPAR solver. With EMPAR, the base case and the contingencies are solved independently, i.e, there is no coupling.

#### Mode
Set SCOPFLOW to either run in `preventive` (0) or `corrective` (1) mode. In preventive mode, any power deficit or surplus in the contingency problem is provided by the swing bus only. In the corrective mode, in addition to the swing bus, thhe generators at PV/PQ buses contribute to the deficit/surplus. The contribution amount is decided by the optimization with the constraint that the real power dispatch for these generators should be within 30-min ramping limit. 

