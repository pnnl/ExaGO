# Unit Test Design for ExaGO OPFLOW's Objective Function

## Goal
Design the scalable unit test for the OPFLOW function **OPFLOWComputeObjective_PBPOL**.

## Objective function

OPFLOW does a minimization of the generation cost and the generation cost function is assumed to be a polynomial function of order 2.
```math
\begin{aligned}
C = \sum_{k=1}^{ng} \alpha_kP^2_{Gk} + \beta_kP_{Gk} + \gamma_k + \sum_{j=1}^{nl}c_{\delta{S_Dj}}({\delta{P^2_{Dj}} + \delta{Q^2_{Dj}}}) 
+ \sum_{i=1}^{nb}c_{\delta{S_i}}({\delta{P^2_{i}} + \delta{Q^2_{i}}})
\end{aligned}
```
where, $`\alpha_k`$, $`\beta_k`$, and $`\gamma_k`$ are the generator $`k`$ cost-cofficients. $`c_{\delta{S_Dj}}`$ is the penalty cost for $`j^{th}`$ load loss, and $`c_{\delta{S_i}}`$ is the penalty cost for power imbalance at bus $`i`$; $`n_g`$, $`n_l`$, and $`n_b`$ are number of generators, number of loads and number of buses; $`P_{Gk}`$ is the real power injection of the generator $`k`$; $`\delta{P_{Dj}}`$ and $`\delta{Q_{Dj}}`$ are real and reactive power load loss for each load $`j`$ if enabled; $`\delta{P^2_{i}}`$ and $`\delta{Q^2_{i}}`$ are real and imaginary power imbalance variables; 

For the purpose of this unit test, last two terms will be assumed zero (and are zero by design).
## Input file
ExaGO OPFLOW reads .m file, thus the input file for unit test is in this format.
A 5-bus system **OF-unittestx1.m** will be used as a basis for this test.

### Parameter values

Following are the value of parameters of interest for this test:

- $`\alpha_{k}=0.01`$
- $`\beta_{k}=0.1`$
- $`\gamma_{k}=8`$
- $`P_{Gk}=10`$

### Solution for N=1

With the parameters of the example network, value of the objective function is equal **10**. 

## Scaling

To build a solution when the network is being multiplied folowing needs to be done:

Multiply 10 with number of segments. (N*10).

### Solution for N=3

For N=3 the bjective function is equal to **30**.
