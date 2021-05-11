# Unit Test Design for ExaGO OPFLOW's  Variable Bounds

## Goal
Design the scalable unit test for the OPFLOW function **OPFLOWSetVariableBounds_PBPOL**.

## Variable bounds
- Bounds on generator k real and reactive power injections: $`P_{Gk}^- \le P_{Gk} \le P_{Gk}^+`$,$`Q_{Gk}^- \le Q_{Gk} \le Q_{Gk}^+`$
- Bounds on real and reactive power load loss for load j: $`0 \le \delta{P_{Dj}} \le P_{Dj}`$,$`0 \le \delta{Q_{Dj}} \le Q_{Dj}`$
- Bounds on voltage magnitude for each bus i: $`(V^-_i) \le V_{i}  \le (V^+_i)`$
- Bound  on the reference bus angle: $`\theta_{ref} = \theta_{ref0}`$, i.e., the reference bus angle is held constant.

## Variable bounds arrays
- 4 bounds per generator
- 4 for load loss
- 2 for voltage per each bus
- 1 for ref angle

ExaGO builds two arrays for the variable bounds (lower and upper). Arrays are built with **for loop** looping through buses and the details of the array elements are given as follows:
<table>
<tr>
<td>Contributor</td> <td>array_l</td> <td>array_u</td>
</tr>
<tr>
<td></td> <td>-INF</td> <td>INF</td>
</tr>
<tr>
<td>bus</td><td>Vmin or 0 (for REF/PV bus)</td> <td>Vmax or 2 (for REF/PV bus)</td>
</tr>
<tr>
<td>gen</td><td>Pgmin</td> <td>Pgmax</td>
</tr>
<tr>
<td>gen</td><td>Qgmin</td> <td>Qgmax</td>
</tr>
</table>


**Determine the order**

## Input file
ExaGO OPFLOW reads .m file, thus the input file for unit test is in this format.
A 5-bus system **SVB-unittestx1.m** will be used as a basis for this test.

### Variable bounds arrays builder for the example network
For the considered network, each array has five bounds for voltages, two bounds for generators and (five more inf bounds???).

<table>
<tr>
<td>array_l</td> <td>-INF</td> <td>Vmin1</td> <td>-INF</td> <td>Vmin2</td> <td>Va*PI/180.0</td> <td>0</td> <td>Pgmin</td> <td>Qgmin</td> <td>-INF</td> <td>Vmin4</td> <td>-INF</td> <td>Vmin5</td>
</tr>
<tr>
<td>array_u</td> <td>INF</td> <td>Vmax1</td> <td>INF</td> <td>Vmax2</td> <td>Va*PI/180.0</td> <td>2</td> <td>Pgmax</td> <td>Qgmax</td> <td>-INF</td> <td>Vmax4</td> <td>-INF</td> <td>Vmax5</td>
</tr>
</table>

### Parameters values
Following are the value of parameters of interest for this test:

- $`V_{min1}=-1`$
- $`V_{min2}=-2`$
- $`V_{min4}=-4`$
- $`V_{min5}=-5`$
- $`V_{max1}=1`$
- $`V_{max2}=2`$
- $`V_{max4}=4`$
- $`V_{max5}=5`$
- $`P_{gmin}=-10`$
- $`Q_{gmin}=-20`$
- $`P_{gmax}=10`$
- $`Q_{gmax}=20`$
- $`V_{a}=57,295779513082320876798154814105`$


### Solution for N=1
With the parameters of the example network arrays are:

<table>
<tr>
<td>array_l</td> <td>-INF</td> <td>-1</td> <td>-INF</td> <td>-2</td> <td>1</td> <td>0</td> <td>-10</td> <td>-20</td> <td>-INF</td> <td>-4</td> <td>-INF</td> <td>-5</td>
</tr>
<tr>
<td>array_u</td> <td>INF</td> <td>1</td> <td>INF</td> <td>2</td> <td>1</td> <td>2</td> <td>10</td> <td>20</td> <td>INF</td> <td>4</td> <td>INF</td> <td>5</td>
</tr>
</table>

## Scaling

To build a solution when the network is being multiplied, to the original array all valuse except first two should be added N times; the elements for the reference ange should be replaced with -INF/INF.

### Solution for N=3
So for N=3 the arrays are:
<table>
<tr>
<td>-INF -1 -INF -2 1 0 -10 -20 -INF -4 -INF -5</td> <td>-INF -2 -INF 0 -10 -20 -INF -4 -INF -5</td> <td>-INF -2 -INF 0 -10 -20 -INF -4 -INF -5</td>
</tr>
<tr>
<td>INF 1 INF 2 1 2 10 20 INF 4 INF 5</td> <td>INF 2 INF 2 10 20 INF 4 INF 5</td> <td>INF 2 INF 2 10 20 INF 4 INF 5</td>
</tr>
</table>
