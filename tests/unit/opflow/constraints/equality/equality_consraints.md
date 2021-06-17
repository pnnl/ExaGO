
# Unit Test Design for ExaGO OPFLOW's  Equality Constraints

## Goal
Design the scalable unit test for the OPFLOW equality constraints.

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

## Input file
### Format
ExaGO OPFLOW reads .m file, thus the input file for unit test is in this format.
### Topology
A network, that consist of:
- 5 buses (1-5)
- 1 generator unit (at bus 3)
- 1 transformer (between buses 2 and 3)
- 3 lines (1-2, 2-4, 4-5)
- 1 switched shunt (at bus 2), and 
- 1 load (at bus 4) is being created.
Figure below shows the oneline diagram of the test network:
![img1.png](one_oneline.jpg)


### Equality constraints array
ExaGO builds array for the equality constraints. Array is build with **for loop** looping through buses and the details of the array elements are given in the following table:
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
<td>2</td> <td>Power imbalance</td> <td>Yes</td> <td>Default is no</td>
</tr>
<tr>
<td>2*num. of gen units</td> <td>Generator</td> <td>No</td> <td></td>
</tr>
<tr>
<td>2*num. of loads points</td> <td>Load</td> <td>No</td> <td>It is calc. diff if the load loss is on</td>
</tr>
<tr>
<td>2*num. of connected lines</td> <td>Line</td> <td>No</td> <td></td>
</tr>
<tr>
<td>2*num. of gen units</td> <td>Generator</td> <td>Yes</td> <td>It is added if gen has a MW setpoint</td>
</tr>
</table>


For the considered network equality constraints vector is as follows:
<table>
<tr>
<td>Contributor</td> <td>Bus1</td> <td>Bus2</td> <td>Bus3</td> <td>Bus4</td> <td>Bus5</td>
</tr>
<tr>
<td>Isolated bus?</td><td>---</td> <td>---</td> <td>---</td> <td>---</td> <td>---</td>
</tr>
<tr>
<td>Shunt</td> <td>0 0</td> <td>$`V_{m2}*V_{m2}*gl`$ $`-V_{m2}*V_{m2}*bl`$</td> <td>0 0</td> <td>0 0</td> <td>0 0</td>
</tr>
<tr>
<td>Power Imbalance</td> <td>---</td> <td>---</td> <td>---</td> <td>---</td> <td>---</td>
</tr>
<tr>
<td>Generators</td> <td>---</td> <td>---</td> <td>$-P_{g}`$ $-Q_{g}`$</td> <td>---</td> <td>---</td>
</tr>
<tr>
<td>Loads</td> <td>---</td> <td>---</td> <td>----</td> <td>$P_{d}`$ $Q_{d}</td> <td>---</td>
</tr>
<tr>
<td>Lines</td> <td>$P_{f12}`$ $Q_{f12}`$</td> <td>$P_{t12}`$ $Q_{t12}`$ $P_{f23}`$ $Q_{f23}`$ $P_{f24}`$ $Q_{f24}`$</td> <td>$P_{t23}`$ $Q_{t23}`$</td> <td>$P_{t24}`$ $Q_{t24}`$ $P_{f45}`$ $Q_{f45}`$</td> <td>$P_{t45}`$ $Q_{t45}`$</td>
</tr>
<tr>
<td>Generators</td> <td>---</td> <td>---</td> <td>---</td> <td>---</td> <td>---</td>
</tr>
</table>



So array has in total 4+8+6+8+4 values.
With the parameters of the example network array is:
<table>
<tr>
<td>Bus1</td> <td>Bus2</td> <td>Bus3</td> <td>Bus4</td> <td>Bus5</td>
</tr>
<tr>
<td>0 0 Pf12 Qf12</td> <td>0 Bl Pt12 Qt12 Pf23 Qf23 Pf24 Qf24</td> <td>0 0 -Pg -Qg Pt34 Qt34</td> <td>0 0 Pd Qd Pt24 Qt24 Pf45 Qf45</td> <td>0 0 Pt45 Qt45 </td>
</tr>
</table>

Where
- $`\P_{f}=G_{ff}*V_{mf}*V_{mf}+V_{mf}*V_{mt}*(G_{ft}*cos(\theta_{ft})+B_{ft}*sin(\theta_{ft}))`$
- $`\Q_{f}=-B_{ff}*V_{mf}*V_{mf}+V_{mf}*V_{mt}*(-B_{ft}*cos(\theta_{ft})+G_{ft}*sin(\theta_{ft}))`$
- $`\P_{t}=G_{tt}*V_{mt}*V_{mt}+V_{mf}*V_{mt}*(G_{tf}*cos(\theta_{tf})+B_{tf}*sin(\theta_{tf}))`$
- $`\Q_{t}=-B_{tt}*V_{mt}*V_{mt}+V_{mf}*V_{mt}*(-B_{ft}*cos(\theta_{tf})+G_{tf}*sin(\theta_{tf}))`$

### Values
Following are the value of parameters of interest for this test:

- $`V_{m}=1`$
- $`\theta=0`$
- $`\R_{branch}=0`$
- $`\X_{branch}=0.000001`$
- $`\B_{branch}=0`$
- $`\tapratio_{branch}=1`$
- $`\phaseshift_{branch}=0`$
- $`P_{g}=10`$
- $`Q_{g}=0`$
- $`P_{d}=10`$
- $`Q_{d}=10`$
- $`G_{l}=0`$
- $`B_{l}=10`$

With this, ** equality constraints array is: 
<table>
<tr>
<td>0 0 -100000 0 0 10 100000 0 -100000 0 -100000 0 0 0 -10 0 100000 0 0 0 10 10 100000 0 -100000 0 0 0 100000 0 </td>
</tr>
</table> **

## Scaling
Idea is to be able to "multiply" the network, and at the same time being able to **easily** evaluate the equality constraints array. Network shown before is considered as a base segment (N=1).

Figure 2 shows the proposed multiplication of the grid with two segments connected (N=2):

![img2.png](two_oneline.jpg)
### Algorithm
Algorithm for the .m file generation with N segments of the base network is:
#### Bus data:
- Copy all but first bus data (column) N times and increment the numbering (First value in the column): (N-1)*4+First.
- Set all Second values in columns of elements with First number euqual to N*4-1 for N>1 (e.g., 7, 11, 15) to 2. This steps ensure that all buses with generator unit are marked as PV buses, and the first one with the generator is SLACK (Secund=3 instead of 2).
- Total number of buses for N segments is N*4+1.
#### Generator data:
- Copy generator data N times.
- ONLY the first value in the generator field (bus number) needs to be updated for each copy = 3+(N-1)*4.
#### Generator cost data:
- Copy generator cost data N times.
#### Branch data:
- Copy all four branches N times.
- First two values in each columns are changing as follows: (N-1)*4+First and (N-1)*4+Second.
#### Bus name data:
- Copy all but first bus name data N times and increment the numbering: (N-1)*4+First. 
#### Generator unit types data:
- Copy generator unit types data N times.
#### Generator fuel types data:
- Copy generator fuel types data N times.

When network is multiplied, following will happen with the array:
Copy once first 4; then for each segment copy the rest of the array and add two more values for all segments except last one.
4+8+6+8+4 +2   +8+6+8+4 +2   +8+6+8+4 +2   +8+6+8+4 +2 +8+6+8+4
Two values added inbetween are Pf and Qf, or -100000 and 0.
So for N=3 the array is:
<table>
<tr>
<td>0 0 -100000 0 0 10 100000 0 -100000 0 -100000 0 0 0 -10 0 100000 0 0 0 10 10 100000 0 -100000 0 0 0 100000 0     -100000 0     0 10 100000 0 -100000 0 -100000 0 0 0 -10 0 100000 0 0 0 10 10 100000 0 -100000 0 0 0 100000 0     -100000 0     0 10 100000 0 -100000 0 -100000 0 0 0 -10 0 100000 0 0 0 10 10 100000 0 -100000 0 0 0 100000 0</td>
</tr>
</table> 