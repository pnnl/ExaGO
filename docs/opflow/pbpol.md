# Optimal Power Flow Formulation in Power Balance Formulation with Polar Representation of Bus Voltages

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

- Two variables at each bus i for the voltage magnitude and angle, $`V_{Mi}`$,$`\theta_{mi}`$. ($`{v}_{mi} = V_{Mi}e^{j\theta_{mi}}`$,$`V_{Mi}= \left | {v}_{mi} \right | `$)
- Two variables at each generator k for the real and reactive power output, $`P_{Gk}`$,$`Q_{Gk}`$
- Two variables at each load j for the real and reactive power load loss, $`\delta{P_{Dj}}`$,$`\delta{Q_{Dj}}`$. These variables are only included 
if the load loss variable flag, -opflow_include_loadloss_variables, is ON.
- Two variables at each bus i for real and imaginary power imbalance (slack) variables, $`\delta{P_{i}}`$,$`\delta{Q_{i}}`$. These variables are included
only if the power imbalance variable flag, -opflow_include_imbalance_variables, is ON.


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

```math
\begin{aligned}
\sum_{A_{br}(f,t) = 1} (G_{ff}(V^2_{Mf}) + V_{Mf}(G_{ft}V_{Mt}\cos(\theta_{mt}) + V_{Mt}B_{ft}\sin(\theta_{mt}))
- \sum_{A_G(f,k) = 1}P_{Gk} + \sum_{A_L(f,j) \neq 0}(P_{Dj} - \delta{P_{Dj}}) + \delta{P_{f}} = \Delta{P_f} = 0 \\
\sum_{A_{br}(f,t) = 1} (-B_{ff}(V^2_{Mf}) + V_{Mf}(G_{ft}V_{Mt}\sin(\theta_{mt})) - V_{Mt}B_{ft}\sin(\theta_{mt}))
    - \sum_{A_G(f,k) \neq 0}Q_{Gk} + \sum_{A_L(f,j) = 1}(Q_{Dj} - \delta{Q_{Dj}}) + \delta{Q_{f}} = \Delta{Q_f} = 0 \\
\end{aligned}
```
Here, $`G_{ff}`$,$`G_{ft}`$ are the self and mutual conductances for line ft, while $`B_{ff}`$,$`B_{ft}`$ are the
self and mutual susceptances, respectively.


## Inequality constraints

### Voltage magnitude constraints for each bus i
The voltage magnitude at bus i is $`V_{Mi} `$. The inequality constraint considered is on the square of the voltage magnitude
```math
(V^-_i)^2 \le V^2_{Mi}  \le (V^+_i)^2
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
## Inequality constraint Jacobian
## Objective Hessian
## Equality constraint Hessian
## Inequality constraint Hessian