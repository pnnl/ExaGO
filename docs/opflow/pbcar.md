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

- Two variables at each bus i for the real and imaginary parts of the complex voltage, $`V_{Ri}`$,$`V_{Ii}`$
- Two variables at each generator k for the real and reactive power output, $`P_{Gk}`$,$`Q_{Gk}`$
- Two variables at each load j for the real and reactive power load loss, $`\Delta{P_{Dj}}`$,$`\Delta{Q_{Dj}}`$. These variables are only included 
if the load loss variable flag, -opflow_include_loadloss_variables, is ON.
- Two variables at each bus i for real and imaginary power imbalance (slack) variables, $`\Delta{P_{i}}`$,$`\Delta{Q_{i}}`$. These variables are included
only if the power imbalance variable flag, -opflow_include_imbalance_variables, is ON.

## Bounds

- Bounds on generator k real and reactive power injections:$`P_{Gk}^- \le P_{Gk} \le P_{Gk}^+`$,$`Q_{Gk}^- \le Q_{Gk} \le Q_{Gk}^+`$
- Bounds on real and reactive power load loss for load j: $`0 \le \Delta{P_{Dj}} \le P_{Dj}`$,$`0 \le \Delta{Q_{Dj}} \le Q_{Dj}`$
- Real and reactive part of complex bus voltage i are unbounded.
- Power imbalance variables are unbounded.


## Objective function

### Minimize generation cost
OPFLOW does a minimization of the generation cost and the generation cost function is assumed to be a polynomial function of order 2.
```math
\begin{aligned}
C = \sum_{k=1}^{ng} \alpha_kP^2_{Gk} + \beta_kP_{Gk} + \gamma_k + \sum_{j=1}^{nl}c_{\Delta{S_Dj}}({\Delta{P^2_{Dj}} + \Delta{Q^2_{Dj}}}) 
+ \sum_{i=1}^{nb}c_{\Delta{S_i}}({\Delta{P^2_{i}} + \Delta{Q^2_{i}}})
\end{aligned}
```
where, $`\alpha_k`$,$`\beta_k`$,$`\gamma_k`$ are the generator $`k`$ cost-cofficients.$`c_{\Delta{S_Dj}}`$ is the penalty cost for jth load loss, 
and $`c_{\Delta{S_i}}`$ is the penalty cost for power imbalance at bus i.

## Equality constraints

### Real and reactive power balance constraints at each bus f
```math
\begin{aligned}
\sum_{A_{br}(f,t) = 1} (G_{ff}(V^2_{Rf} + V^2_{If}) + V_{Rf}(G_{ft}V_{Rt} - B_{ft}V_{It}) + V_{If}(B_{ft}V_{Rt} + G_{ft}V_{It}))
    - \sum_{A_G(f,k) = 1}P_{Gk} + \sum_{A_L(f,j) \neq 0}(P_{Dj} - \Delta{P_{Dj}}) + \Delta{P_{f}} = 0 \\
\sum_{A_{br}(f,t) = 1} (-B_{ff}(V^2_{Rf} + V^2_{If}) + V_{If}(G_{ft}V_{Rt} - B_{ft}V_{It}) - V_{Rf}(B_{ft}V_{Rt} + G_{ft}V_{It}))
    - \sum_{A_G(f,k) \neq 0}Q_{Gk} + \sum_{A_L(f,j) = 1}(Q_{Dj} - \Delta{Q_{Dj}}) + \Delta{Q_{f}} = 0 \\
\end{aligned}
```
### Voltage angle constraint at ref. bus
This constraint attempts to hold the reference bus angle fixed at $`\angle{\bar{V_{ref}}}`$.
```math
\begin{aligned}
    V_{Iref} - V_{Rref}\text{tan}(\angle{\bar{V_{ref}}}) = 0
\end{aligned}
```

## Inequality constraints

### Voltage magnitude constraints for each bus i
```math
(V^-_i)^2 \le V^2_{Ri} + V^2_{Ii} \le (V^+_i)^2
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
0 \le S_f \le (S^+_{ft})^2 \\
0 \le S_t \le (S^+_{ft})^2
\end{aligned}
```
where the maximum flow,$`S^+_{ft}`$ is either the RATE_A (normal), RATE_B (short-term), or RATE_C (emergency) rating of the line.

## Gradient

```math
\begin{aligned}
\dfrac{\partial{C}}{\partial{P_{Gk}}} &= 2\alpha_kP_{Gk} + \beta_k \\
\dfrac{\partial{C}}{\partial{\Delta{P_{Dj}}}} &= 2c_{\Delta{S_{Dj}}}\Delta{P_{Dj}} \\
\dfrac{\partial{C}}{\partial{\Delta{Q_{Dj}}}} &= 2c_{\Delta{S_{Dj}}}\Delta{Q_{Dj}} \\
\dfrac{\partial{C}}{\partial{\Delta{P_{i}}}} &= 2c_{\Delta{S_i}}\Delta{P_{i}} \\
\dfrac{\partial{C}}{\partial{\Delta{Q_{i}}}} &= 2c_{\Delta{S_i}}\Delta{Q_{i}}

\end{aligned}
```


## Equality constraint Jacobian
## Inequality constraint Jacobian
## Objective Hessian
## Equality constraint Hessian
## Inequality constraint Hessian