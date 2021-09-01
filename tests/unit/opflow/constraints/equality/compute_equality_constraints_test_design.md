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
self and mutual susceptances, respectively and are calculated form the line parameters as follows:

- $`G_{ff}=y_{ff}[0]=\dfrac{R}{R^2+X^2}\dfrac{1}{tap^2}`$
- $`B_{ff}=y_{ff}[1]=(\dfrac{-X}{R^2+X^2}+\dfrac{B_{c}}{2})\dfrac{1}{tap^2}`$
- $`G_{ft}=y_{ft}[0]=(-\dfrac{R}{R^2+X^2}*tap*\cos(shift)-\dfrac{X}{R^2+X^2}*tap*\sin(shift))\dfrac{1}{tap^2}`$
- $`B_{ft}=y_{ft}[1]=(\dfrac{X}{R^2+X^2}*tap*\cos(shift)-\dfrac{R}{R^2+X^2}*tap*\sin(shift))\dfrac{1}{tap^2}`$
- $`G_{tf}=y_{tf}[0]=(-\dfrac{R}{R^2+X^2}*tap*\cos(shift)+\dfrac{X}{R^2+X^2}*tap*\sin(shift))\dfrac{1}{tap^2}`$
- $`B_{tf}=y_{tf}[1]=(\dfrac{X}{R^2+X^2}*tap*\cos(shift)+\dfrac{R}{R^2+X^2}*tap*\sin(shift))\dfrac{1}{tap^2}`$
- $`G_{tt}=y_{tt}[0]=\dfrac{R}{R^2+X^2}`$
- $`B_{tt}=y_{tt}[1]=\dfrac{-X}{R^2+X^2}+\dfrac{B_{c}}{2}`$

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
<td>2 * num. of gen units</td> <td>Generator</td> <td>Yes</td> <td>Values will be added if the option *Generator real power setpoint set? ** is ON. Default is OFF</td>
</tr>
</table>

## Input

ExaGO OPFLOW reads .m file, thus the input file for unit test is in this format.
A 5-bus system **CEC-unittestx1.m** will be used as a basis for this test. In addition, an artifical solution vector will be also generated as an input for the test.

### Parameters values in .m file

Following are the values of parameters of interest for this test:
- $`R_{branch}=2`$
- $`X_{branch}=1`$
- $`B_{branch}=1.2`$
- $`tapratio_{branch}=1`$ for all but transformer that has $`tapratio_{transformer}=2`$
- $`phaseshift_{branch}=0`$ for all but transformer that has $`phaseshift_{branch}=60`$
- $`P_{d}=-3.4`$
- $`Q_{d}=8.8`$
- $`G_{l}=0.25`$
- $`B_{l}=-0.05`$

### Solution vector values

Following are the values of the solution vector's elements of interest for this test:

- $`theta=0`$ for all but transformer that has $`theta=30`$
- $`V_{m}=2`$
- $`P_{g}=1.6`$
- $`Q_{g}=-2.2`$

## Vectors for N=1

### Solution vector for N=1

Solution vector can be built with **OPFLOWSetInitialGuess_PBPOL**, but for the purpose of this unit test, solution vector will be build independently in such a way that only the values of the interest will be different than zero. 
In general, solution vector has following elements per bus:
1. Voltage angle
2. Voltage magnitude
3. Generator MW (if generator bus)
4. Generator MVar (if generator bus)

For the 5-bus system **CEC-unittestx1.m**, solution vector is:
<table>
<tr>
<td>0</td> <td>2</td> <td>0</td> <td>2</td> <td>30</td> <td>2</td> <td>1.6</td> <td>-2.2</td> <td>0</td> <td>2</td> <td>0</td> <td>2</td> 
</tr>
</table>

### Residual vector builder for equality constraints for the example network

For the considered network residual vector for equality constraints is built as follows:
<table>
<tr>
<td>Contributor</td> <td>Bus1[0]</td> <td>Bus1[1]</td> <td>Bus2[0]</td> <td>Bus2[1]</td> <td>Bus3[0]</td> <td>Bus3[1]</td> <td>Bus4[0]</td> <td>Bus4[1]</td> <td>Bus5[0]</td> <td>Bus5[1]</td>
</tr>
<tr>
<td>Isolated bus?</td><td>-</td> <td>-</td> <td>-</td> <td>-</td> <td>-</td> <td>-</td> <td>-</td> <td>-</td> <td>-</td> <td>-</td>
</tr>
<tr>
<td>Shunt</td> <td>0</td> <td>0</td> <td> Vm2*Vm2*Gl</td> <td>-Vm2*Vm2*Bl</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td>
</tr>
<tr>
<td>Power Imbalance</td> <td>-</td> <td>-</td> <td>-</td> <td>-</td> <td>-</td> <td>-</td> <td>-</td> <td>-</td> <td>-</td> <td>-</td> 
</tr>
<tr>
<td>Generators</td> <td>-</td> <td>-</td>  <td>-</td> <td>-</td> <td>-Pg</td>  <td>-Qg</td> <td>-</td> <td>-</td> <td>-</td> <td>-</td>
</tr>
<tr>
<td>Loads</td> <td>-</td> <td>-</td> <td>-</td> <td>-</td> <td>-</td> <td>-</td> <td>Pd</td> <td>Qd</td> <td>-</td> <td>-</td>
</tr>
<tr>
<td>Lines</td> <td>Pf12</td> <td>Qf12</td> <td>Pt12+Pf23+Pf24</td> <td>Qt12+Qf23+Qf24</td> <td>Pt23</td> <td>Qt23</td> <td>Pt24+Pf45</td> <td>Qt24+ Qf45</td> <td>Pt45</td> <td>Qt45</td>
</tr>
<tr>
<td>Generators</td> <td>-</td> <td>-</td> <td>-</td> <td>-</td> <td>-</td> <td>-</td> <td>-</td> <td>-</td> <td>-</td> <td>-</td>
</tr>
</table>

So array has in total 2+2+2+2+2 values.

<table>
<tr>
<td>Bus1[0]</td> <td>Bus1[1]</td> <td>Bus2[0]</td> <td>Bus2[1]</td> <td>Bus3[0]</td> <td>Bus3[1]</td> <td>Bus4[0]</td> <td>Bus4[1]</td> <td>Bus5[0]</td> <td>Bus5[1]</td>
</tr>
<tr>
<td>Pf12</td> <td>Qf12</td> <td>Vm2*Vm2*Gl+Pt12+Pf23+Pf24</td> <td>-Vm2*Vm2*Bl+Qt12+Qf23+Qf24</td> <td>-Pg+Pt23</td> <td>-Qg+Qt23</td> <td>Pd+Pt24+Pf45</td> <td>Qd+Qt24+Qf45</td> <td>Pt45</td> <td>Qt45</td>
</tr>
</table>

Where
- $`P_{f}=G_{ff}*V_{mf}*V_{mf}+V_{mf}*V_{mt}*(G_{ft}*cos(\theta_{ft})+B_{ft}*sin(\theta_{ft}))`$
- $`Q_{f}=-B_{ff}*V_{mf}*V_{mf}+V_{mf}*V_{mt}*(-B_{ft}*cos(\theta_{ft})+G_{ft}*sin(\theta_{ft}))`$
- $`P_{t}=G_{tt}*V_{mt}*V_{mt}+V_{mf}*V_{mt}*(G_{tf}*cos(\theta_{tf})+B_{tf}*sin(\theta_{tf}))`$
- $`Q_{t}=-B_{tt}*V_{mt}*V_{mt}+V_{mf}*V_{mt}*(-B_{ft}*cos(\theta_{tf})+G_{tf}*sin(\theta_{tf}))`$



### Residual vector for the equality constraints for N=1

With the parameters of the example network the vector is:

<table>
<tr>
<td>2.8</td> <td>0.0</td> <td>8.8</td> <td>2.2</td> <td>2.2</td> <td>2.2</td> <td>2.2</td> <td>8.8</td> <td>2.8</td> <td>0.0</td>
</tr>
</table>

## Scaling

### Solution vector

To scale the solution vector, following needs to be done:
1. Copy once original solution vector.
2. Copy N-1 times all elements of the initial vector, except first two.

### Solution vector for N=3

For N=3 the solution vector is:
<table>
<tr>
<td>0   2   0   2   30   2   1.6   -2.2   0   2   0   2</td> <td>0   2   30   2   1.6   -2.2   0   2   0   2</td><td>0   2   30   2   1.6   -2.2   0   2   0   2</td>
</tr>
</table>

### Residual vector

To build a residual vector when the network is being multiplied folowing needs to be done:

Copy only once first 2 elements; next 6 elements will be copied without any modifications for each segment; last two elements are first multipled by 2 and then coppied for each N except last one where you copy without multiplication.

Example:
a b c d e f g h i j

 a b only once

 c d e f g h

 ix2 jx2

 c d e f g h

 ix2 jx2 

 c d e f g h

 i j  

### Residual for N=3

So for N=3 the vector is:
<table>
<tr>
<td>2.8  0.0</td> <td> 8.8  2.2  2.2  2.2  2.2  8.8</td> <td> 5.6  0.0</td><td> 8.8  2.2  2.2  2.2  2.2  8.8</td> <td> 5.6  0.0</td><td> 8.8  2.2  2.2  2.2  2.2  8.8</td> <td> 2.8  0.0</td>
</tr>
</table>
