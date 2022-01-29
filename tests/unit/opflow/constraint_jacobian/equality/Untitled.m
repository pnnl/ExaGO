Vm=2;
Vm1=2;
Vm2=2;
Vm3=2;
Vm4=2;
Vm5=2;

theta1=deg2rad(0);
theta2=deg2rad(0);
theta3=deg2rad(30);
theta4=deg2rad(0);
theta5=deg2rad(0);


R=2;
X=1;
Bc=1.2;

tap=1;
tapt=2;
shift=deg2rad(0);
shiftt=deg2rad(60);


Pg=1.6;
Qg=-2.2;

Pd=-3.4;
Qd=8.8;

Gl=0.25;
Bl=-0.05;


theta12=theta1-theta2;
theta21=theta2-theta1;

theta23=theta2-theta3;
theta32=theta3-theta2;

theta24=theta2-theta4;
theta42=theta4-theta2;

theta45=theta4-theta5;
theta54=theta5-theta4;

Zm=R*R+X*X;
G=R/Zm;
B=-X/Zm;

tap2=tap*tap;
tapt2=tapt*tapt;

tapr = tap*cos(shift);
tapi = tap*sin(shift);
taprt = tapt*cos(shiftt);
tapit = tapt*sin(shiftt);

Gff = G/tap2; 
Bff = (B+Bc/2.0)/tap2;        
Gft = -(G*tapr - B*tapi)/tap2;
Bft = -(B*tapr + G*tapi)/tap2;
Gtf = -(G*tapr + B*tapi)/tap2;
Btf = -(B*tapr - G*tapi)/tap2;
Gtt = G;
Btt = B+Bc/2.0;

Gfft = G/tapt2; 
Bfft = (B+Bc/2.0)/tapt2;        
Gftt = -(G*taprt - B*tapit)/tapt2;
Bftt = -(B*taprt + G*tapit)/tapt2;
Gtft = -(G*taprt + B*tapit)/tapt2;
Btft = -(B*taprt - G*tapit)/tapt2;
Gttt = G;
Bttt = B+Bc/2.0; 



Pf12 = Gff*Vm1*Vm1+Vm1*Vm2*( Gft*cos(theta12)+Bft*sin(theta12));
Qf12 =-Bff*Vm1*Vm1+Vm1*Vm2*(-Bft*cos(theta12)+Gft*sin(theta12));
Pt12 = Gtt*Vm2*Vm2+Vm1*Vm2*( Gtf*cos(theta21)+Btf*sin(theta21));
Qt12 =-Btt*Vm2*Vm2+Vm1*Vm2*(-Btf*cos(theta21)+Gtf*sin(theta21));

dPf12dVm1=2*Gff*Vm1+Vm2*(Gft*cos(theta12)+Bft*sin(theta12));
dPf12dVm2=Vm2*(Gft*cos(theta12)+Bft*sin(theta12));
dPf12dtheta1=Vm1*Vm2*(-Gft*sin(theta12)+Bft*cos(theta12));
dPf12dtheta2=Vm1*Vm2*(Gft*sin(theta12)-Bft*cos(theta12));

dQf12dVm1=-2*Bff*Vm1+Vm2*(-Bft*cos(theta12)+Gft*sin(theta12));
dQf12dVm2=Vm1*(-Bft*cos(theta12)+Gft*sin(theta12));
dQf12dtheta1=Vm1*Vm2*(Bft*sin(theta12)+Gft*cos(theta12));
dQf12dtheta2=Vm1*Vm2*(-Bft*sin(theta12)-Gft*cos(theta12));

dPt12dVm1 = Vm2*( Gtf*cos(theta21)+Btf*sin(theta21));
dPt12dVm2 = 2*Gtt*Vm2+Vm1*( Gtf*cos(theta21)+Btf*sin(theta21));
dPt12dtheta1 = Vm1*Vm2*( Gtf*sin(theta21)-Btf*cos(theta21));
dPt12dtheta2 = Vm1*Vm2*(-Gtf*sin(theta21)+Btf*cos(theta21));

dQt12dVm1 =Vm2*(-Btf*cos(theta21)+Gtf*sin(theta21));
dQt12dVm2 =-2*Btt*Vm2+Vm1*(-Btf*cos(theta21)+Gtf*sin(theta21));
dQt12dtheta1 =Vm1*Vm2*(-Btf*sin(theta21)-Gtf*cos(theta21));
dQt12dtheta2 =Vm1*Vm2*(Btf*sin(theta21)+Gtf*cos(theta21));


Pf23 = Gfft*Vm2*Vm2+Vm2*Vm3*( Gftt*cos(theta23)+Bftt*sin(theta23));
Qf23 =-Bfft*Vm2*Vm2+Vm2*Vm3*(-Bftt*cos(theta23)+Gftt*sin(theta23));
Pt23 = Gttt*Vm3*Vm3+Vm2*Vm3*( Gtft*cos(theta32)+Btft*sin(theta32));
Qt23 =-Bttt*Vm3*Vm3+Vm2*Vm3*(-Btft*cos(theta32)+Gtft*sin(theta32));

dPf23dVm2=2*Gfft*Vm2+Vm3*(Gftt*cos(theta23)+Bftt*sin(theta23));
dPf23dVm3=Vm3*(Gftt*cos(theta23)+Bftt*sin(theta23));
dPf23dtheta2=Vm2*Vm3*(-Gftt*sin(theta23)+Bftt*cos(theta23));
dPf23dtheta3=Vm2*Vm3*(Gftt*sin(theta23)-Bftt*cos(theta23));

dQf23dVm2=-2*Bfft*Vm2+Vm3*(-Bftt*cos(theta23)+Gftt*sin(theta23));
dQf23dVm3=Vm2*(-Bftt*cos(theta23)+Gftt*sin(theta23));
dQf23dtheta2=Vm2*Vm3*(Bftt*sin(theta23)+Gftt*cos(theta23));
dQf23dtheta3=Vm2*Vm3*(-Bftt*sin(theta23)-Gftt*cos(theta23));

dPt23dVm2 = Vm3*( Gtft*cos(theta32)+Btft*sin(theta32));
dPt23dVm3 = 2*Gttt*Vm3+Vm2*( Gtft*cos(theta32)+Btft*sin(theta32));
dPt23dtheta2 = Vm2*Vm3*( Gtft*sin(theta32)-Btft*cos(theta32));
dPt23dtheta3 = Vm2*Vm3*(-Gtft*sin(theta32)+Btft*cos(theta32));

dQt23dVm2 =Vm3*(-Btft*cos(theta32)+Gtft*sin(theta32));
dQt23dVm3 =-2*Bttt*Vm3+Vm2*(-Btft*cos(theta32)+Gtft*sin(theta32));
dQt23dtheta2 =Vm2*Vm3*(-Btft*sin(theta32)-Gtft*cos(theta32));
dQt23dtheta3 =Vm2*Vm3*(Btft*sin(theta32)+Gtft*cos(theta32));

Pf24 = Gff*Vm2*Vm2+Vm2*Vm4*( Gft*cos(theta24)+Bft*sin(theta24));
Qf24 =-Bff*Vm2*Vm2+Vm2*Vm4*(-Bft*cos(theta24)+Gft*sin(theta24));
Pt24 = Gtt*Vm4*Vm4+Vm2*Vm4*( Gtf*cos(theta42)+Btf*sin(theta42));
Qt24 =-Btt*Vm4*Vm4+Vm2*Vm4*(-Btf*cos(theta42)+Gtf*sin(theta42));

dPf24dVm2=2*Gff*Vm2+Vm4*(Gft*cos(theta24)+Bft*sin(theta24));
dPf24dVm4=Vm4*(Gft*cos(theta24)+Bft*sin(theta24));
dPf24dtheta2=Vm2*Vm4*(-Gft*sin(theta24)+Bft*cos(theta24));
dPf24dtheta4=Vm2*Vm4*(Gft*sin(theta24)-Bft*cos(theta24));

dQf24dVm2=-2*Bff*Vm2+Vm4*(-Bft*cos(theta24)+Gft*sin(theta24));
dQf24dVm4=Vm2*(-Bft*cos(theta24)+Gft*sin(theta24));
dQf24dtheta2=Vm2*Vm4*(Bft*sin(theta24)+Gft*cos(theta24));
dQf24dtheta4=Vm2*Vm4*(-Bft*sin(theta24)-Gft*cos(theta24));

dPt24dVm2 = Vm4*( Gtf*cos(theta42)+Btf*sin(theta42));
dPt24dVm4 = 2*Gtt*Vm4+Vm2*( Gtf*cos(theta42)+Btf*sin(theta42));
dPt24dtheta2 = Vm2*Vm4*( Gtf*sin(theta42)-Btf*cos(theta42));
dPt24dthera4 = Vm2*Vm4*(-Gtf*sin(theta42)+Btf*cos(theta42));

dQt24dVm2 =Vm4*(-Btf*cos(theta42)+Gtf*sin(theta42));
dQt24dVm4 =-2*Btt*Vm4+Vm2*(-Btf*cos(theta42)+Gtf*sin(theta42));
dQt24dtheta2 =Vm2*Vm4*(-Btf*sin(theta42)-Gtf*cos(theta42));
dQt24dtheta4 =Vm2*Vm4*(Btf*sin(theta42)+Gtf*cos(theta42));

Pf45 = Gff*Vm4*Vm4+Vm4*Vm5*( Gft*cos(theta45)+Bft*sin(theta45));
Qf45 =-Bff*Vm4*Vm4+Vm4*Vm5*(-Bft*cos(theta45)+Gft*sin(theta45));
Pt45 = Gtt*Vm5*Vm5+Vm4*Vm5*( Gtf*cos(theta54)+Btf*sin(theta54));
Qt45 =-Btt*Vm5*Vm5+Vm4*Vm5*(-Btf*cos(theta54)+Gtf*sin(theta54));

dPf45dVm4=2*Gff*Vm4+Vm5*(Gft*cos(theta45)+Bft*sin(theta45));
dPf45dVm5=Vm5*(Gft*cos(theta45)+Bft*sin(theta45));
dPf45dtheta4=Vm4*Vm5*(-Gft*sin(theta45)+Bft*cos(theta45));
dPf45dtheta5=Vm4*Vm5*(Gft*sin(theta45)-Bft*cos(theta45));

dQf45dVm4=-2*Bff*Vm4+Vm5*(-Bft*cos(theta45)+Gft*sin(theta45));
dQf45dVm5=Vm4*(-Bft*cos(theta45)+Gft*sin(theta45));
dQf45dtheta4=Vm4*Vm5*(Bft*sin(theta45)+Gft*cos(theta45));
dQf45dtheta5=Vm4*Vm5*(-Bft*sin(theta45)-Gft*cos(theta45));

dPt45dVm4 = Vm5*( Gtf*cos(theta54)+Btf*sin(theta54));
dPt45dVm5 = 2*Gtt*Vm5+Vm4*( Gtf*cos(theta54)+Btf*sin(theta54));
dPt45dtheta4 = Vm4*Vm5*( Gtf*sin(theta54)-Btf*cos(theta54));
dPt45dtheta5 = Vm4*Vm5*(-Gtf*sin(theta54)+Btf*cos(theta54));

dQt45dVm4 =Vm5*(-Btf*cos(theta54)+Gtf*sin(theta54));
dQt45dVm5 =-2*Btt*Vm5+Vm4*(-Btf*cos(theta54)+Gtf*sin(theta54));
dQt45dtheta4 =Vm4*Vm5*(-Btf*sin(theta54)-Gtf*cos(theta54));
dQt45dtheta5 =Vm4*Vm5*(Btf*sin(theta54)+Gtf*cos(theta54));

% number of buses
N = 5;
% number of generators
Gen = 1;
J=zeros(2*N, 2*N+2*Gen)
J(1,1) = dPf12dtheta1;
J(2,1) = dQf12dtheta1;
J(3,1) = dPt12dtheta1;
J(4,1) = dQt12dtheta1;

J(1,2) = dPf12dVm1;
J(2,2) = dQf12dVm1;
J(3,2) = dPt12dVm1;
J(4,2) = dQt12dVm1;

J(1,3) = dPf12dtheta2;
J(2,3) = dQf12dtheta2;
J(3,3) = dPt12dtheta2+0+dPf23dtheta2+dPf24dtheta2;
J(4,3) = dQt12dtheta2+0+dQf23dtheta2+dQf24dtheta2;
J(5,3) = dPt23dtheta2;
J(6,3) = dQt23dtheta2;
J(7,3) = dPt24dtheta2;
J(8,3) = dQt24dtheta2;

J(1,4) = dPf12dVm2;
J(2,4) = dQf12dVm2;
J(3,4) = dPt12dVm2+2*Vm*Gl+dPf23dVm2+dPf24dVm2;
J(4,4) = dQt12dVm2-2*Vm*Bl+dQf23dVm2+dQf24dVm2;
J(5,4) = dPt23dVm2;
J(6,4) = dQt23dVm2;
J(7,4) = dPt24dVm2;
J(8,4) = dQt24dVm2;

J(3,5) = dPf23dtheta3;
J(4,5) = dQf23dtheta3;
J(5,5) = dPt23dtheta3;
J(6,5) = dQt23dtheta3;

J(3,6) = dPf23dVm3;
J(4,6) = dQf23dVm3;
J(5,6) = dPt23dVm3;
J(6,6) = dQt23dVm3;

J(5,7) = -1;
J(6,8) = -1;

J(3,9) = dPf24dtheta4;
J(4,9) = dQf24dtheta4;
J(7,9) = dPt24dthera4+dPf45dtheta4;
J(8,9) = dQt24dtheta4+dQf45dtheta4;
J(9,9) = dPt45dtheta4;
J(10,9) = dQt45dtheta4;

J(3,10) = dPf24dVm4;
J(4,10) = dQf24dVm4;
J(7,10) = dPt24dVm4+dPf45dVm4;
J(8,10) = dQt24dVm4+dQf45dVm4;
J(9,10) = dPt45dVm4;
J(10,10) = dQt45dVm4;

J(7,11) = dPf45dtheta5;
J(8,11) = dQf45dtheta5;
J(9,11) = dPt45dtheta5;
J(10,11) = dQt45dtheta5;

J(7,12) = dPf45dVm5;
J(8,12) = dQf45dVm5;
J(9,12) = dPt45dVm5;
J(10,12) = dQt45dVm5;

J












