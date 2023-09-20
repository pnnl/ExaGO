# Unit Test Design for ExaGO OPFLOW's  Equality Constraint Jacobian

## Goal
Design the scalable unit test for the OPFLOW function **OPFLOWComputeEqualityConstraintJacobian_PBPOL**.

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

## Equality constraint Jacobian

Jacobian terms are calculated as partial w.r.t voltage magniude and angle at each bus, as well as w.r.t generator injections, (and load loss and power imbalance if enabled).
Dimensions of the matrix are:
- $`Number_of_rows = 2*Number_of_buses`
- $`Number_of_columns = 2*(Number_of_buses+Number_of_generators)+(2*Number_of_buses if load loss and power imbalance options are enabled)`

## Input

\ExaGO OPFLOW reads .m file, thus the input file for unit test is in this format.
A 5-bus system **CECJ-unittestx1.m** will be used as a basis for this test. In addition, an artifical solution vector will be also generated as an input for the test.

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

For the 5-bus system **CECJ-unittestx1.m**, solution vector is:
<table>
<tr>
<td>0</td> <td>2</td> <td>0</td> <td>2</td> <td>30</td> <td>2</td> <td>1.6</td> <td>-2.2</td> <td>0</td> <td>2</td> <td>0</td> <td>2</td> 
</tr>
</table>

### Matrix builder for equality constraints Jacobian for the example network

For the considered network matrix for equality constraint Jacobian is built as follows:
<table>
<tr>
<td> </td> <td>theta1</td> <td>Vm1</td> <td>theta2</td> <td>Vm2</td> <td>theta3</td> <td>Vm3</td> <td>Pg3</td> <td>Qg3</td> <td>theta4</td> <td>Vm4</td> <td>theta5</td> <td>Vm5</td>
</tr>
<tr>
<td>P1</td> <td>dPf12dtheta1</td> <td>dPf12dVm1</td> <td>dPf12dtheta2</td> <td>dPf12dVm2</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td>
</tr>
<tr>
<td>Q1</td> <td>dQf12dtheta1</td> <td>dQf12dVm1</td> <td>dQf12dtheta2</td> <td>dQf12dVm2</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td>
</tr>
<tr>
<td>P2</td> <td>dPt12dtheta1</td> <td>dPt12dVm1</td> <td>dPt12dtheta2+0+dPf23dtheta2+dPf24dtheta2</td> <td>dPt12dVm2+2*Vm*Gl+dPf23dVm2+dPf24dVm2</td> <td>dPf23dtheta3</td> <td>dPf23dVm3</td> <td>0</td> <td>0</td> <td>dPf24dtheta4</td> <td>dPf24dVm4</td> <td>0</td> <td>0</td>
</tr>
<tr>
<td>Q2</td> <td>dQt12dtheta1</td> <td>dQt12dVm1</td> <td>dQt12dtheta2+0+dQf23dtheta2+dQf24dtheta2</td> <td>dQt12dVm2+2*Vm*Bl+dQf23dVm2+dQf24dVm2</td> <td>dQf23dtheta3</td> <td>dQf23dVm3</td> <td>0</td> <td>0</td> <td>dQf24dtheta4</td> <td>dQf24dVm4</td> <td>0</td> <td>0</td>
</tr>
<tr>
<td>P3</td> <td>0</td> <td>0</td> <td>dPt23dtheta2</td> <td>dPt23dVm2</td> <td>dPt23dtheta3</td> <td>dPt23dVm3</td> <td>-1</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td>
</tr>
<tr>
<td>Q3</td> <td>0</td> <td>0</td> <td>dQt23dtheta2</td> <td>dQt23dVm2</td> <td>dQt23dtheta3</td> <td>dQt23dVm3</td> <td>0</td> <td>-1</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td>
</tr>
<tr>
<td>P4</td> <td>0</td> <td>0</td> <td>dPt24dtheta2</td> <td>dPt24dVm2</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>dPt24dthera4+dPf45dtheta4</td> <td>dPt24dVm4+dPf45dVm4</td> <td>dPf45dtheta5</td> <td>dPf45dVm5</td>
</tr>
<tr>
<td>Q4</td> <td>0</td> <td>0</td> <td>dQt24dtheta2</td> <td>dQt24dVm2</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>dQt24dthera4+dQf45dtheta4</td> <td>dQt24dVm4+dQf45dVm4</td> <td>dQf45dtheta5</td> <td>dQf45dVm5</td>
</tr>
<tr>
<td>P5</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>dPt45dtheta4</td> <td>dPt45dVm4</td> <td>dPt45dtheta5</td> <td>dPt45dVm5</td>
</tr>
<tr>
<td>Q5</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>dQt45dtheta4</td> <td>dQt45dVm4</td> <td>dQt45dtheta5</td> <td>dQt45dVm5</td>
</tr>
</table>


### Jacobian matrix for N=1

With the parameters of the example network the matrix is:

<table>
<tr>
<td> </td> <td>theta1</td> <td>Vm1</td> <td>theta2</td> <td>Vm2</td> <td>theta3</td> <td>Vm3</td> <td>Pg3</td> <td>Qg3</td> <td>theta4</td> <td>Vm4</td> <td>theta5</td> <td>Vm5</td>
</tr>
<tr>
<td>P1</td> <td>0.8</td> <td>0.8</td> <td>-0.8</td> <td>-0.8</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td>
</tr>
<tr>
<td>Q1</td> <td>-1.6</td> <td>-2.0</td> <td>1.6</td> <td>-0.4</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td>
</tr>
<tr>
<td>P2</td> <td>-0.8</td> <td>-0.8</td> <td>0.8</td> <td>2.8</td> <td>0.8</td> <td>-0.2</td> <td>0</td> <td>0</td> <td>-0.8</td> <td>-0.8</td> <td>0</td> <td>0</td>
</tr>
<tr>
<td>Q2</td> <td>1.6</td> <td>-0.4</td> <td>-3.6</td> <td>-3.8</td> <td>0.4</td> <td>0.4</td> <td>0</td> <td>0</td> <td>1.6</td> <td>-0.4</td> <td>0</td> <td>0</td>
</tr>
<tr>
<td>P3</td> <td>0</td> <td>0</td> <td>-0.8</td> <td>0.2</td> <td>0.8</td> <td>1.8</td> <td>-1</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td>
</tr>
<tr>
<td>Q3</td> <td>0</td> <td>0</td> <td>-0.4</td> <td>-0.4</td> <td>0.4</td> <td>-2.0</td> <td>0</td> <td>-1</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td>
</tr>
<tr>
<td>P4</td> <td>0</td> <td>0</td> <td>-0.8</td> <td>-0.8</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>1.6</td> <td>1.6</td> <td>-0.8</td> <td>-0.8</td>
</tr>
<tr>
<td>Q4</td> <td>0</td> <td>0</td> <td>1.6</td> <td>-0.4</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>-3.2</td> <td>-4.0</td> <td>1.6</td> <td>-0.4</td>
</tr>
<tr>
<td>P5</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>-0.8</td> <td>-0.8</td> <td>0.8</td> <td>0.8</td>
</tr>
<tr>
<td>Q5</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>1.6</td> <td>-0.4</td> <td>-1.6</td> <td>-2.0</td>
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

### Jacobian

To build a Jacobian matrix when the network is being multiplied the process for matrix builder is shown in the Figure below:

![img1.png](Jacobian.jpg)

### Jacobian for N=3
