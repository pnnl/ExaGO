function mpc = case9mod
%CASE9a    Power flow data for 9 bus, 3 generator case.
%   Please see CASEFORMAT for details on the case file format.
%
% Based on data from Sauer and Pai's book 'Power System Dynamics and
% Stability' Chapter 7
%
% Modifications
% Gen. 3 is modeled as a wind generator
% Changed bus 5 reactive power load to 30 MVAr

%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100.0000;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
    1   3    0.0   0.0  0.0 0.000000 1   1.0400  0.0000  345.0000    1   1.02  0.98;
    2   2    0.0   0.0  0.0 0.000000 1   1.0250  9.3000  345.0000    1   1.02  0.98;
    3   2    0.0   0.0  0.0 0.000000 1   1.0250  4.7000  345.0000    1   1.02  0.98;
    4   1    0.0   0.0  0.0      0.0 1   1.0260 -2.2000  345.0000    1   1.1000  0.9000;
    5   1   75.0  30.0  0.0 0.000000 1   0.9960 -4.0000  345.0000    1   1.1000  0.9000;
    6   1   90.0  30.0  0.0 0.000000 1   1.0130 -3.7000  345.0000    1   1.1000  0.9000;
    7   1    0.0   0.0  0.0      0.0 1   1.0260  3.7000  345.0000    1   1.1000  0.9000;
    8   1  100.0  35.0  0.0 0.000000 1   1.0160  0.7000  345.0000    1   1.1000  0.9000;
    9   1    0.0   0.0  0.0      0.0 1   1.0320  2.0000  345.0000    1   1.1000  0.9000;
];

%% generator data
    %	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
    1   71.6000  27.0000  300.0000  -300.0000  1.0000  100.0000  1   350.0000    10.0000 0	0	0	0	0	0	0	0	0	0	0;
    2  163.0000   6.7000  300.0000  -300.0000  1.00  100.0000  1   300.0000    10.0000 0	0	0	0	0	0	0	0	0	0	0;
    3   45.0000 -10.9000  0.0000     0.0000  1.00  100.0000  1   75.0000      75.0000 0	0	0	0	0	0	0	0	0	0	0;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
    1   4   0.0001  0.0576  0.0001  380.0000    250.0000    250.0000    0.0000  0.0000  1 -360	360;
    2   7   0.0001  0.0625  0.0001  250.0000    250.0000    250.0000    0.0000  0.0000  1 -360	360;
    3   9   0.0001  0.0586  0.0001  300.0000    300.0000    250.0000    0.0000  0.0000  1 -360	360;
    4   5   0.0100  0.0850  0.1760  250.0000    250.0000    250.0000    0.0000  0.0000  1 -360	360;
    4   6   0.0170  0.0920  0.1580  250.0000    250.0000    250.0000    0.0000  0.0000  1 -360	360;
    5   7   0.0320  0.1610  0.3060  250.0000    250.0000    250.0000    0.0000  0.0000  1 -360	360;
    6   9   0.0390  0.1700  0.3580  150.0000    150.0000    150.0000    0.0000  0.0000  1 -360	360;
    7   8   0.0085  0.0720  0.1490  250.0000    250.0000    250.0000    0.0000  0.0000  1 -360	360;
    8   9   0.0119  0.1008  0.2090  150.0000    150.0000    150.0000    0.0000  0.0000  1 -360	360;
];

%%-----  OPF Data  -----%%

%% area data
%	area	refbus
mpc.areas = [
	1	1;
];

%% generator cost data
%	1	startup	shutdown	n	x1	y1	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
	2	1500.00	0.00	3	0.11	5	150;
	2	2000.00	0.00	3	0.085	1.2	600;
	2	0000.00	0.00	3       0.0     0.0      0.0;
];

%% generator unit type (see GENTYPES)
mpc.gentype = {
	'ST';
	'ST';
	'W2';
};

%% generator fuel type (see GENFUELS)
mpc.genfuel = {
	'coal';
	'coal';
        'wind';
};


return;
