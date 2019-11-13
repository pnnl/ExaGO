# Optimal Power Flow Formulation in Current Balance Formulation with Cartesian Representation of Bus Voltages

## Sizes
<table>
<tr>
<td>Number of buses</td> <td>n_b</td>
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

## Bounds

- Bounds on generator k real and reactive power injections:$`P_{Gk}^- \le P_{Gk} \le P_{Gk}^+`$,$`Q_{Gk}^- \le Q_{Gk} \le Q_{Gk}^+`$
- Bounds on real and reactive power load loss for load j: $`0 \le \delta{P_{Dj}} \le P_{Dj}`$,$`0 \le \delta{Q_{Dj}} \le Q_{Dj}`$



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