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

## Input file
ExaGO OPFLOW reads .m file, thus the input file for unit test is in this format.
A 5-bus system **OFG-unittestx1.m** will be used as a basis for this test.

### Parameter values

Following are the value of parameters of interest for this test:

- $`\alpha_{k}=0.045`$
- $`\beta_{k}=0.1`$
- $`P_{Gk}=10`$

### Solution for N=1

With the parameters of the example network, value of the objective function is equal **1**. 

## Scaling

To build a solution when the network is being multiplied folowing needs to be done:

Create the vector with lenght(vector)=N.

### Solution for N=3

For N=3 the value of the objective function's gradient is equal to vector with length of 3 where all values are equal to 1:

<table>
<tr>
<td>1 1 1 </td>
</tr>
</table> 

