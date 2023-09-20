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

## Input

ExaGO OPFLOW reads .m file, thus the input file for unit test is in this format.
A 5-bus system **OF-unittestx1.m** will be used as a basis for this test. In addition, an artifical solution vector will be also generated as an input for the test.

### Parameter values in .m file

Following are the values of parameters of interest for this test:

- $`\alpha_{k}=0.01`$
- $`\beta_{k}=0.1`$
- $`\gamma_{k}=8`$

### Solution vector values

Following is the value of the solution vector's element of interest for this test:

- $`P_{Gk}=10`$

## Vectors for N=1

### Solution vector for N=1

Solution vector can be built with **OPFLOWSetInitialGuess_PBPOL**, but for the purpose of this unit test, solution vector will be build independently in such a way that only the values of the interest will be different than zero. 
In general, solution vector has following elements per bus:
1. Voltage angle
2. Voltage magnitude
3. Generator MW (if generator bus)
4. Generator MVar (if generator bus)

For the 5-bus system **OF-unittestx1.m**, solution vector is:
<table>
<tr>
<td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>10</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> 
</tr>
</table>

### Residual for the objective function for N=1

With the parameters of the example network, value of the residual of the objective function is equal **10**. 

## Scaling

### Solution vector

To scale the solution vector, following needs to be done:
1. Copy once original solution vector.
2. Copy N-1 times all elements of the initial vector, except first two.

### Solution vector for N=3

For N=3 the solution vector is:
<table>
<tr>
<td>0   0   0   0   0   0   10   0   0   0   0   0</td> <td>0   0   0   0   10   0   0   0   0   0</td><td>0   0   0   0   10   0   0   0   0   0</td>
</tr>
</table>

### Residual vector

To calculate a residual when the network is being multiplied following needs to be done:
1. Multiply 10 with the number of segments. (N*10).

### Residual for N=3

For N=3 the residual of objective function is equal to **30**.
