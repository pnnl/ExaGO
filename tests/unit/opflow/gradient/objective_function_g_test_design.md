# Unit Test Design for ExaGO OPFLOW's Objective Function's Gradient

## Goal
Design the scalable unit test for the OPFLOW function **OPFLOWComputeGradient_PBPOL**.

## Objective function's gradient

OPFLOW does a minimization of the generation cost and the generation cost function $`C`$ is assumed to be a polynomial function of order 2 which gradient is calculated as:
```math
\begin{aligned}
\dfrac{\partial{C}}{\partial{P_{Gk}}} &= 2\alpha_kP_{Gk} + \beta_k \\
\dfrac{\partial{C}}{\partial{\delta{P_{Dj}}}} &= 2c_{\delta{S_{Dj}}}\delta{P_{Dj}} \\
\dfrac{\partial{C}}{\partial{\delta{Q_{Dj}}}} &= 2c_{\delta{S_{Dj}}}\delta{Q_{Dj}} \\
\dfrac{\partial{C}}{\partial{\delta{P_{i}}}} &= 2c_{\delta{S_i}}\delta{P_{i}} \\
\dfrac{\partial{C}}{\partial{\delta{Q_{i}}}} &= 2c_{\delta{S_i}}\delta{Q_{i}}
\end{aligned}
```
where, $`\alpha_k`$, and $`\beta_k`$, are the generator $`k`$ cost-cofficients. $`c_{\delta{S_Dj}}`$ is the penalty cost for $`j^{th}`$ load loss, and $`c_{\delta{S_i}}`$ is the penalty cost for power imbalance at bus $`i`$; $`P_{Gk}`$ is the real power injection of the generator $`k`$; $`\delta{P_{Dj}}`$ and $`\delta{Q_{Dj}}`$ are real and reactive power load loss for each load $`j`$ if enabled; $`\delta{P_{i}}`$ and $`\delta{Q_{i}}`$ are real and imaginary power imbalance variables; 

For the purpose of this unit test, focus is only on the first gradient, because the other four are not calculated by default and are also zero if included. Those elements are include if following flags are enabled:
- include_loadloss_variables
- include_powerimbalance_variables

## Input

ExaGO OPFLOW reads .m file, thus the input file for unit test is in this format.
A 5-bus system **OFG-unittestx1.m** will be used as a basis for this test. In addition, an artifical solution vector will be also generated as an input for the test.

### Parameter values in .m file

Following are the values of parameters of interest for this test:

- $`\alpha_{k}=0.045`$
- $`\beta_{k}=0.1`$

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

For the 5-bus system **OFG-unittestx1.m**, solution vector is:
<table>
<tr>
<td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>10</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> 
</tr>
</table>

### Residual for the gradient for N=1

With the parameters of the example network, value of the gradient is equal **1**. 

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

To build a residual vector when the network is being multiplied folowing needs to be done:

1. Create the vector with lenght(vector)=N.

### Residual vector for N=3

For N=3 residual vector of gradient is equal to vector with length of 3 where all values are equal to 1:

<table>
<tr>
<td>1 1 1 </td>
</tr>
</table> 

