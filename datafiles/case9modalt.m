function mpc = case9modalt
%CASE9MODALT

%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	3	0	0	0	0	1	1	0	345	1	1.1	0.9;
	2	1	0	0	0	0	1	0.957565778	0.756823816	345	1	1.1	0.9;
	3	2	0	0	0	0	1	1	21.3936369	345	1	1.1	0.9;
	4	1	0	0	0	0	1	0.973103168	-1.3656722	345	1	1.1	0.9;
	5	1	100	50	0	0	1	0.939690211	-4.03130362	345	1	1.1	0.9;
	6	1	100	30	0	0	1	0.949511987	-0.299521257	345	1	1.1	0.9;
	7	1	0	0	0	0	1	0.957562785	0.756824102	345	1	1.1	0.9;
	8	1	100	35	0	0	1	0.949187274	2.99825809	345	1	1.1	0.9;
	9	1	0	0	0	0	1	0.983558715	12.1397484	345	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
	1	40.3459846	47.1007254	300	-300	1	100	1	350	10	0	0	0	0	0	0	0	0	0	0	0;
	2	0	           0	    300	-300	1	100	0	300	10	0	0	0	0	0	0	0	0	0	0	0;
	3	269.991593	49.4350388	300	-300	1	100	1	270	10	0	0	0	0	0	0	0	0	0	0	0;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax	Pf	Qf	Pt	Qt
mpc.branch = [
	1	4	0.0001	0.0576	0.0001	380	250	250	0	0	1	-360	360	40.3460	47.1007	-40.3421	-44.8947;
	2	7	0.0001	0.0625	0.0001	250	250	250	0	0	1	-360	360	0.0000	0.0000	-0.0000	-0.0092;
	3	9	0.0001	0.0586	0.0001	300	300	250	0	0	1	-360	360	269.9916	49.4350	-269.9163	-5.2958;
	4	5	0.01	0.085	0.176	250	250	250	0	0	1	-360	360	53.9226	24.7393	-53.5000	-37.2510;
	4	6	0.017	0.092	0.158	250	250	250	0	0	1	-360	360	-13.5804	20.1554	13.7507	-33.8374;
	5	7	0.032	0.161	0.306	250	250	250	0	0	1	-360	360	-46.5000	-12.7490	47.2838	-10.8467;
	6	9	0.039	0.17	0.358	150	150	150	0	0	1	-360	360	-113.7507	3.8374	119.5205	-12.1413;
	7	8	0.0085	0.072	0.149	250	250	250	0	0	1	-360	360	-47.2838	10.8559	47.5200	-22.3979;
	8	9	0.0119	0.1008	0.209	150	150	150	0	0	1	-360	360	-147.5200	-12.6021	150.3958	17.4371;
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
	2	1500	0	3	100	4	150;
	2	1500	0	3	100	4	150;
	2	1500	0	3	100	4	150;
];
