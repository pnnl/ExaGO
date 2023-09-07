## Multi-period optimal power flow(TCOPFLOW)
TCOPFLOW solves a multi-period AC optimal power flow problem with the objective of minimizing the total cost over the given time horizon while adhering to intra-time and inter-temporal constraints. The problem setup is as following:

```math
\begin{aligned}
\text{min}&~\sum_{t \in N_t} f(x_t)& \\
&\text{s.t.}& \\
&~g(x_t) = 0~~~t \in \{0,N_t\}& \\
&~h(x_t) \le 0~~t \in \{0,N_t\}& \\
&x^- \le x_t \le x^+~~t\in \{0,N_t\}& \\
-\Delta{x}_t & \le x_t - x_{t-\Delta{t}} \le \Delta{x}_t~~t \in \{1,N_t\}&
\end{aligned}
```

Here, $`N_t`$ is the number of timesteps. Each time-step has a full AC optimal power flow formulation with equality constraints $`g(x_t)`$, inequality $`h(x_t)`$ and the lower/upper bounds $`x^-, x^+`$. Additionally, the temporal coupling constraints, given by the last equation, need to be satisfied. The temporal coupling constraints in TCOPFLOW are the generator real power output ramping limits between successive time-steps.

### Dependency
ExaGO must be built with Ipopt for this application.

### Usage
TCOPFLOW is executed via
```
mpiexec -n <N> ./tcopflow <options>
```

where \<options\> are the available command line options as given in the next section.

### Options
There are several options available for the current version of TCOPFLOW. These options can be either set through the options file `options/tcopflowoptions` or via the command line.

#### Network file (-netfile \<netfilename\>): 
Set the name of the network file. Only MATPOWER format is currently supported. 4096 characters max.

```
mpiexec -n <N> ./tcopflow -netfile <netfilename>
```

#### Wind generation profile (-tcopflow_windgenprofile \<windgenprofile_filename\>)
The name of the file describing the wind generation profile. See `$EXAGO_DIR/datafiles/case9/case9mod_gen3_wind.m` for the format.

#### Active power load profile (-tcopflow_ploadprofile \<ploadprofile_filename\>)
The name of the file describing the active load profile. See `$EXAGO_DIR/datafiles/case9/load_P.csv` for the format.

#### Reactive power load profile (-tcopflow_qloadprofile \<qloadprofile_filename\>)
The name of the file describing the reactive load profile. See `$EXAGO_DIR/datafiles/case9/load_Q.csv` for the format.

#### Time-step (-tcopflow_dT \<time_step\>)
The time-step for multi-period TCOPFFLOW in minutes.

#### Duration (-tcopflow_duration \<duration\>)
The duration for the multi-period TCOFLOW in hours.
