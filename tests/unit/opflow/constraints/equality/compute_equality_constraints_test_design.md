# Unit Test Design for ExaGO OPFLOW's  Equality Constraints

## Goal
Design the scalable unit test for the OPFLOW function **OPFLOWComputeEqualityConstraints_PBPOL**.

## Equality constraints

```math
\begin{aligned}
\sum_{A_{br}(f,t) = 1} (G_{ff}(V^2_{f}) + V_{f}V_{t}(G_{ft}\cos(\theta_{f}-\theta_t) + B_{ft}\sin(\theta_{f}-\theta_t))
- \sum_{A_G(f,k) = 1}P_{Gk} + \sum_{A_L(f,j) \neq 0}(P_{Dj} - \delta{P_{Dj}}) + \delta{P_{f}} = \Delta{P_f} = 0 \\
\sum_{A_{br}(f,t) = 1} (-B_{ff}(V^2_{f}) + V_{f}V_{t}(G_{ft}\sin(\theta_{f}-\theta_t) - B_{ft}\cos(\theta_{f}-\theta_t))
    - \sum_{A_G(f,k) \neq 0}Q_{Gk} + \sum_{A_L(f,j) = 1}(Q_{Dj} - \delta{Q_{Dj}}) + \delta{Q_{f}} = \Delta{Q_f} = 0 \\
\end{aligned}
```
Here, $`G_{ff}`$,$`G_{ft}`$ are the self and mutual conductances for line ft, while $`B_{ff}`$,$`B_{ft}`$ are the
self and mutual susceptances, respectively.

## Equality constraints array
ExaGO builds array of the equality constraints. Array is build with **for loop** looping through buses and the details of the array elements are given in the following table:
<table>
<tr>
<td>Number of values</td> <td>Contributor</td> <td>Optional?</td> <td>Comment</td>
</tr>
<tr>
<td>2</td> <td>Isolated bus?</td> <td>Yes</td> <td>Values will be added only if the bus is isolated</td>
</tr>
<tr>
<td>2</td> <td>Shunt</td> <td>No</td> <td></td>
</tr>
<tr>
<td>2</td> <td>Power imbalance</td> <td>Yes</td> <td>Values will be added if the option **Allow power imbalance at buses** is ON. Default is OFF</td>
</tr>
<tr>
<td>2 * num. of gen units</td> <td>Generator</td> <td>No</td> <td></td>
</tr>
<tr>
<td>2 * num. of loads points</td> <td>Load</td> <td>No</td> <td>Values will be diffrently calculated if the option **Include load loss** is ON. Default is OFF</td>
</tr>
<tr>
<td>2 * num. of connected lines</td> <td>Line</td> <td>No</td> <td></td>
</tr>
<tr>
<td>2 * num. of gen units</td> <td>Generator</td> <td>Yes</td> <td>Values will be added if the option *Generator real power set point set? ** is ON. Default is OFF</td>
</tr>
</table>

## Input file
ExaGO OPFLOW reads .m file, thus the input file for unit test is in this format.
A 5-bus system **CEC-unittestx1.m** will be used as a basis for this test.

### Equality constraints array builder for the example network

For the considered network equality constraints vector is built as follows:
<table>
<tr>
<td>Contributor</td> <td>Bus1</td> <td>Bus2</td> <td>Bus3</td> <td>Bus4</td> <td>Bus5</td>
</tr>
<tr>
<td>Isolated bus?</td><td>---</td> <td>---</td> <td>---</td> <td>---</td> <td>---</td>
</tr>
<tr>
<td>Shunt</td> <td>0 0</td> <td>Vm2*Vm2*Gl -Vm2*Vm2*Bl</td> <td>0 0</td> <td>0 0</td> <td>0 0</td>
</tr>
<tr>
<td>Power Imbalance</td> <td>---</td> <td>---</td> <td>---</td> <td>---</td> <td>---</td>
</tr>
<tr>
<td>Generators</td> <td>---</td> <td>---</td> <td>-Pg -Qg</td> <td>---</td> <td>---</td>
</tr>
<tr>
<td>Loads</td> <td>---</td> <td>---</td> <td>----</td> <td>Pd Qd</td> <td>---</td>
</tr>
<tr>
<td>Lines</td> <td>Pf12 Qf12</td> <td>Pt12 Qt12 Pf23 Qf23 Pf24 Qf24</td> <td>Pt23 Qt23</td> <td>Pt24 Qt24 Pf45 Qf45</td> <td>Pt45 Qt45</td>
</tr>
<tr>
<td>Generators</td> <td>---</td> <td>---</td> <td>---</td> <td>---</td> <td>---</td>
</tr>
</table>

So array has in total 4+8+6+8+4 values.

<table>
<tr>
<td>Bus1</td> <td>Bus2</td> <td>Bus3</td> <td>Bus4</td> <td>Bus5</td>
</tr>
<tr>
<td>0 0 Pf12 Qf12</td> <td>0 -Bl Pt12 Qt12 Pf23 Qf23 Pf24 Qf24</td> <td>0 0 -Pg -Qg Pt34 Qt34</td> <td>0 0 Pd Qd Pt24 Qt24 Pf45 Qf45</td> <td>0 0 Pt45 Qt45 </td>
</tr>
</table>

Where
- $`P_{f}=G_{ff}*V_{mf}*V_{mf}+V_{mf}*V_{mt}*(G_{ft}*cos(\theta_{ft})+B_{ft}*sin(\theta_{ft}))`$
- $`Q_{f}=-B_{ff}*V_{mf}*V_{mf}+V_{mf}*V_{mt}*(-B_{ft}*cos(\theta_{ft})+G_{ft}*sin(\theta_{ft}))`$
- $`P_{t}=G_{tt}*V_{mt}*V_{mt}+V_{mf}*V_{mt}*(G_{tf}*cos(\theta_{tf})+B_{tf}*sin(\theta_{tf}))`$
- $`Q_{t}=-B_{tt}*V_{mt}*V_{mt}+V_{mf}*V_{mt}*(-B_{ft}*cos(\theta_{tf})+G_{tf}*sin(\theta_{tf}))`$

### Parameters values
Following are the value of parameters of interest for this test:

- $`V_{m}=1`$
- $`theta=0`$
- $`R_{branch}=0`$
- $`X_{branch}=0.000001`$
- $`B_{branch}=0`$
- $`tapratio_{branch}=1`$
- $`phaseshift_{branch}=0`$
- $`P_{g}=10`$
- $`Q_{g}=0`$
- $`P_{d}=10`$
- $`Q_{d}=10`$
- $`G_{l}=0`$
- $`B_{l}=10`$

### Solution for N=1
With the parameters of the example network array is:

<table>
<tr>
<td>0 0 -100000 0 0 -10 100000 0 -100000 0 -100000 0 0 0 -10 0 100000 0 0 0 10 10 100000 0 -100000 0 0 0 100000 0 </td>
</tr>
</table>

## Scaling

To build a solution when the network is being multiplied folowing needs to be done:

Copy once first 4 elements; then for each segment (N) copy the rest of the array and add two more values for all segments except last one.
Two values added inbetween are Pf and Qf, or -100000 and 0.


4+8+6+8+4 +2   +8+6+8+4 +2   +8+6+8+4

### Solution for N=3
So for N=3 the array is:
<table>
<tr>
<td>0 0 -100000 0 0 -10 100000 0 -100000 0 -100000 0 0 0 -10 0 100000 0 0 0 10 10 100000 0 -100000 0 0 0 100000 0     -100000 0     0 -10 100000 0 -100000 0 -100000 0 0 0 -10 0 100000 0 0 0 10 10 100000 0 -100000 0 0 0 100000 0     -100000 0     0 -10 100000 0 -100000 0 -100000 0 0 0 -10 0 100000 0 0 0 10 10 100000 0 -100000 0 0 0 100000 0</td>
</tr>
</table> 
