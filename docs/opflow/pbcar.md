# Optimal Power Flow Formulation in Power Balance Formulation with Cartesian Representation of Bus Voltages

## Sizes
<table>
<tr>
<td>Number of buses</td> <td>nb</td>
</tr>
<tr>
<td>Number of branches</td> <td>nbr</td>
</tr>
<tr>
<td>Number of generators</td> <td>ng</td>
</tr>
<tr>
<td>Number of loads</td> <td>nl</td>
</tr>
</table>

## Variables

- Two variables at each bus i for the real and imaginary parts of the complex voltage, $`V_{Ri}`$,$`V_{Ii}`$. ($`\bar{V}_i = V_{Ri} + \sqrt{-1}V_{Ii}`$)
- Two variables at each generator k for the real and reactive power output, $`P_{Gk}`$,$`Q_{Gk}`$
- Two variables at each load j for the real and reactive power load loss, $`\delta{P_{Dj}}`$,$`\delta{Q_{Dj}}`$. These variables are only included 
if the load loss variable flag, -opflow_include_loadloss_variables, is ON.
- Two variables at each bus i for real and imaginary power imbalance (slack) variables, $`\delta{P_{i}}`$,$`\delta{Q_{i}}`$. These variables are included
only if the power imbalance variable flag, -opflow_include_imbalance_variables, is ON.

## Bounds

- Bounds on generator k real and reactive power injections:$`P_{Gk}^- \le P_{Gk} \le P_{Gk}^+`$,$`Q_{Gk}^- \le Q_{Gk} \le Q_{Gk}^+`$
- Bounds on real and reactive power load loss for load j: $`0 \le \delta{P_{Dj}} \le P_{Dj}`$,$`0 \le \delta{Q_{Dj}} \le Q_{Dj}`$
- Real and reactive part of complex bus voltage i are unbounded.
- Power imbalance variables are unbounded.


## Objective function

### Minimize generation cost
OPFLOW does a minimization of the generation cost and the generation cost function is assumed to be a polynomial function of order 2.
```math
\begin{aligned}
C = \sum_{k=1}^{ng} \alpha_kP^2_{Gk} + \beta_kP_{Gk} + \gamma_k + \sum_{j=1}^{nl}c_{\delta{S_Dj}}({\delta{P^2_{Dj}} + \delta{Q^2_{Dj}}}) 
+ \sum_{i=1}^{nb}c_{\delta{S_i}}({\delta{P^2_{i}} + \delta{Q^2_{i}}})
\end{aligned}
```
where, $`\alpha_k`$,$`\beta_k`$,$`\gamma_k`$ are the generator $`k`$ cost-cofficients.$`c_{\delta{S_Dj}}`$ is the penalty cost for jth load loss, 
and $`c_{\delta{S_i}}`$ is the penalty cost for power imbalance at bus i.

## Equality constraints

### Real and reactive power balance constraints at each bus f
```math
\begin{aligned}
\sum_{A_{br}(f,t) = 1} (G_{ff}(V^2_{Rf} + V^2_{If}) + V_{Rf}(G_{ft}V_{Rt} - B_{ft}V_{It}) + V_{If}(B_{ft}V_{Rt} + G_{ft}V_{It}))
    - \sum_{A_G(f,k) = 1}P_{Gk} + \sum_{A_L(f,j) \neq 0}(P_{Dj} - \delta{P_{Dj}}) + \delta{P_{f}} = \Delta{P_f} = 0 \\
\sum_{A_{br}(f,t) = 1} (-B_{ff}(V^2_{Rf} + V^2_{If}) + V_{If}(G_{ft}V_{Rt} - B_{ft}V_{It}) - V_{Rf}(B_{ft}V_{Rt} + G_{ft}V_{It}))
    - \sum_{A_G(f,k) \neq 0}Q_{Gk} + \sum_{A_L(f,j) = 1}(Q_{Dj} - \delta{Q_{Dj}}) + \delta{Q_{f}} = \Delta{Q_f} = 0 \\
\end{aligned}
```
Here, $`G_{ff}`$,$`G_{ft}`$ are the self and mutual conductances for line ft, while $`B_{ff}`$,$`B_{ft}`$ are the
self and mutual susceptances, respectively.
### Voltage angle constraint at ref. bus
This constraint attempts to hold the reference bus angle fixed at $`\theta_{ref}`$.
```math
\begin{aligned}
    V_{Iref} - V_{Rref}\text{tan}(\theta_{ref}) = \Delta{\theta_{ref}} = 0
\end{aligned}
```

## Inequality constraints

### Voltage magnitude constraints for each bus i
The voltage magnitude at bus i is $`V_i = \sqrt{V^2_{Ri} + V^2_{Ii}}`$. The inequality constraint considered is on the square of the voltage magnitude
```math
(V^-_i)^2 \le V^2_i = V^2_{Ri} + V^2_{Ii} \le (V^+_i)^2
```

### Apparent line flow constraints for each line between 'from' bus f and 'to' bus t
Let $`P_{ft}`$ and $`Q_{ft}`$ be the real and reactive power flow from bus f to bus t on line ft. 
Similarly, let $`P_{tf}`$ and $`Q_{tf}`$ be the real and reactive power flow from bus t to bus f.
The apparent power flows $`S_{f}`$ and $`S_{t}`$ at the from and to ends of the line are given by
```math
\begin{aligned}
S_f = \sqrt{P^2_{ft} + Q^2_{ft}} \\
S_t = \sqrt{P^2_{tf} + Q^2_{tf}}
\end{aligned}
```
The apparent line flow constraints are considered on the square of the from and to bus flows.
```math
\begin{aligned}
0 \le S^2_f \le (S^+_{ft})^2 \\
0 \le S^2_t \le (S^+_{ft})^2
\end{aligned}
```
where the maximum flow,$`S^+_{ft}`$ is either the RATE_A (normal), RATE_B (short-term), or RATE_C (emergency) rating of the line.

## Gradient

```math
\begin{aligned}
\dfrac{\partial{C}}{\partial{P_{Gk}}} &= 2\alpha_kP_{Gk} + \beta_k \\
\dfrac{\partial{C}}{\partial{\delta{P_{Dj}}}} &= 2c_{\delta{S_{Dj}}}\delta{P_{Dj}} \\
\dfrac{\partial{C}}{\partial{\delta{Q_{Dj}}}} &= 2c_{\delta{S_{Dj}}}\delta{Q_{Dj}} \\
\dfrac{\partial{C}}{\partial{\delta{P_{i}}}} &= 2c_{\delta{S_i}}\delta{P_{i}} \\
\dfrac{\partial{C}}{\partial{\delta{Q_{i}}}} &= 2c_{\delta{S_i}}\delta{Q_{i}}
\end{aligned}
```

## Equality constraint Jacobian

### Jacobian terms for power balance constraints at bus f

#### Partials with respect to from and to bus voltage variables
```math
\begin{aligned}
\dfrac{\partial{\Delta{P_f}}}{\partial{V_{Rf}}} &= \sum_{A_{br}(f,t) = 1} 2G_{ff}V_{Rf} + G_{ft}V_{Rt} - B_{ft}V_{it}\\
\dfrac{\partial{\Delta{P_f}}}{\partial{V_{If}}} &= \sum_{A_{br}(f,t) = 1} 2G_{ff}V_{If} + B_{ft}V_{Rt} + G_{ft}V_{it}\\
\dfrac{\partial{\Delta{P_f}}}{\partial{V_{Rt}}} &= \sum_{A_{br}(f,t) = 1}                 G_{ft}V_{Rf} + B_{ft}V_{If}\\
\dfrac{\partial{\Delta{P_f}}}{\partial{V_{It}}} &= \sum_{A_{br}(f,t) = 1}                -B_{ft}V_{Rf} + G_{ft}V_{If}
\end{aligned}
```

```math
\begin{aligned}
\dfrac{\partial{\Delta{Q_f}}}{\partial{V_{Rf}}} &= \sum_{A_{br}(f,t) = 1} -2B_{ff}V_{Rf} - B_{ft}V_{Rt} - G_{ft}V_{it}\\
\dfrac{\partial{\Delta{Q_f}}}{\partial{V_{If}}} &= \sum_{A_{br}(f,t) = 1} -2B_{ff}V_{If} + G_{ft}V_{Rt} - B_{ft}V_{it}\\
\dfrac{\partial{\Delta{Q_f}}}{\partial{V_{Rt}}} &= \sum_{A_{br}(f,t) = 1}                  G_{ft}V_{If} - B_{ft}V_{Rf}\\
\dfrac{\partial{\Delta{Q_f}}}{\partial{V_{It}}} &= \sum_{A_{br}(f,t) = 1}                 -B_{ft}V_{If} - G_{ft}V_{Rf}
\end{aligned}
```

#### Partials w.r.t generator k injections, load loss for load j, and power imbalance at bus f

```math
\begin{aligned}
\dfrac{\partial{\Delta{P_f}}}{\partial{P_{Gk}}} &= -1 \\
\dfrac{\partial{\Delta{Q_f}}}{\partial{Q_{Gk}}} &= -1 \\
\dfrac{\partial{\Delta{P_f}}}{\partial{\delta{P_{Dj}}}} &= -1 \\
\dfrac{\partial{\Delta{Q_f}}}{\partial{\delta{Q_{Dj}}}} &= -1 \\
\dfrac{\partial{\Delta{P_f}}}{\partial{\delta{P_{i}}}} &= 1 \\
\dfrac{\partial{\Delta{Q_f}}}{\partial{\delta{Q_{i}}}} &= 1
\end{aligned}
```

### Jacobian elements for voltage angle constraint at ref. bus

### Voltage angle constraint at ref. bus

```math
\begin{aligned}
    \dfrac{\partial{\Delta{\theta_{ref}}}}{\partial{V_{Rref}}} &= -\text{tan}(\theta_{ref}) \\
    \dfrac{\partial{\Delta{\theta_{ref}}}}{\partial{V_{Iref}}} &= 1
\end{aligned}
```
## Inequality constraint Jacobian

### Jacobian elements for voltage magniitude constraint at bus i
```math
\begin{aligned}
\dfrac{\partial{V^2_i}}{\partial{V_{Ri}}} = 2V_{Ri} \\
\dfrac{\partial{V^2_i}}{\partial{V_{Ii}}} = 2V_{Ii}
\end{aligned}
```

### Jacobian elements for line ft flow inequality constraint
The from and to bus real and reactive power flows on line ft are

```math
\begin{aligned}
  P_f &=  G_{ff}(V^2_{Rf} + V^2_{If}) + V_{Rf}(G_{ft}V_{Rt} - B_{ft}V_{It}) + V_{If}(B_{ft}V_{Rt} + G_{ft}V_{It}) \\
  Q_f &= -B_{ff}(V^2_{Rf} + V^2_{If}) + V_{If}(G_{ft}V_{Rt} - B_{ft}V_{It}) - V_{Rf}(B_{ft}V_{Rt} + G_{ft}V_{It}) \\
  P_t &=  G_{tt}(V^2_{Rt} + V^2_{It}) + V_{Rt}(G_{tf}V_{Rf} - B_{tf}V_{If}) + V_{It}(B_{tf}V_{Rf} + G_{tf}V_{If}) \\
  Q_t &= -B_{tt}(V^2_{Rt} + V^2_{It}) + V_{It}(G_{tf}V_{Rf} - B_{tf}V_{If}) - V_{Rt}(B_{tf}V_{Rf} + G_{tf}V_{If})
\end{aligned}
```

```math
\begin{aligned}
\dfrac{\partial{S^2_f}}{\partial{V_{Rf}}} &= \dfrac{\partial{S^2_f}}{\partial{P_f}}\dfrac{\partial{P_f}}{\partial{V_{Rf}}} 
                                           + \dfrac{\partial{S^2_f}}{\partial{Q_f}}\dfrac{\partial{Q_f}}{\partial{V_{Rf}}}
                                          &= 2P_f(2G_{ff}V_{Rf} + G_{ft}V_{Rt} - B_{ft}V_{it}) + 2Q_f(-2B_{ff}V_{Rf} - B_{ft}V_{Rt} - G_{ft}V_{it})
\end{aligned}
```

Similarly,
```math
\begin{aligned}
\dfrac{\partial{S^2_f}}{\partial{V_{If}}} &= 2P_f(2G_{ff}V_{If} + B_{ft}V_{Rt} + G_{ft}V_{it}) + 2Q_f(-2B_{ff}V_{If} + G_{ft}V_{Rt} - B_{ft}V_{it}) \\
\dfrac{\partial{S^2_f}}{\partial{V_{Rt}}} &= 2P_f(G_{ft}V_{Rf} + B_{ft}V_{If}) + 2Q_f(G_{ft}V_{If} - B_{ft}V_{Rf})\\
\dfrac{\partial{S^2_f}}{\partial{V_{It}}} &= 2P_f(-B_{ft}V_{Rf} + G_{ft}V_{If}) + 2Q_f(-B_{ft}V_{If} - G_{ft}V_{Rf})
\end{aligned}
```

```math
\begin{aligned}
\dfrac{\partial{S^2_t}}{\partial{V_{Rt}}} &= \dfrac{\partial{S^2_t}}{\partial{P_t}}\dfrac{\partial{P_t}}{\partial{V_{Rt}}} 
                                           + \dfrac{\partial{S^2_t}}{\partial{Q_t}}\dfrac{\partial{Q_t}}{\partial{V_{Rt}}}
                                          &= 2P_t(2G_{tt}V_{Rt} + G_{tf}V_{Rf} - B_{tf}V_{if}) + 2Q_t(-2B_{tt}V_{Rt} - B_{tf}V_{Rf} - G_{tf}V_{if})
\end{aligned}
```

Similarly,
```math
\begin{aligned}
\dfrac{\partial{S^2_t}}{\partial{V_{It}}} &= 2P_t(2G_{tt}V_{It} + B_{tf}V_{Rf} + G_{tf}V_{if}) + 2Q_t(-2B_{tt}V_{It} + G_{tf}V_{Rf} - B_{tf}V_{if}) \\
\dfrac{\partial{S^2_t}}{\partial{V_{Rf}}} &= 2P_t(G_{tf}V_{Rt} + B_{tf}V_{It}) + 2Q_t(G_{tf}V_{It} - B_{tf}V_{Rt})\\
\dfrac{\partial{S^2_t}}{\partial{V_{If}}} &= 2P_t(-B_{tf}V_{Rt} + G_{tf}V_{It}) + 2Q_t(-B_{tf}V_{It} - G_{tf}V_{Rt})
\end{aligned}
```



## Objective Hessian
## Equality constraint Hessian
## Inequality constraint Hessian