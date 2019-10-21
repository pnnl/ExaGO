# Optimal Power Flow Formulation in Power Balance Formulation with Polar Representation of Bus Voltages

## Variables
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

## Gradient
## Equality constraints
## Inequality constraints
## Equality constraint Jacobian
## Inequality constraint Jacobian
## Objective Hessian
## Equality constraint Hessian
## Inequality constraint Hessian