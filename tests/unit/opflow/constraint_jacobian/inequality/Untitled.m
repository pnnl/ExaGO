Vm=2;
Vm1=2;
Vm2=2;
Vm3=2;
Vm4=2;
Vm5=2;

Vset = 0.5;
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
Qmax=197.8;
Qmin=-202.2;

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

IEC1 = (Qg - Qmax)*(Vset-Vm);
IEC2 = (Qmin - Qg)*(Vm-Vset);

dIEC1dQg = Vset-Vm;
dIEC2dQg = Vset-Vm;

dIEC1dVm3 = Qmax-Qg;
dIEC2dVm3 = Qmin-Qg;

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

Sf12 = Pf12*Pf12 + Qf12*Qf12;
St12 = Pt12*Pt12 + Qt12*Qt12;

dSf12dPf12 = 2*Pf12;
dSf12dQf12 = 2*Qf12;
dSt12dPt12 = 2*Pt12;
dSt12dQt12 = 2*Qt12;

multiplier = 1;

dSf12dVm1 =    (dSf12dPf12*dPf12dVm1+   dSf12dQf12*dQf12dVm1)*multiplier;
dSf12dVm2 =    (dSf12dPf12*dPf12dVm2+   dSf12dQf12*dQf12dVm2)*multiplier;
dSf12dtheta1 = (dSf12dPf12*dPf12dtheta1+dSf12dQf12*dQf12dtheta1)*multiplier;
dSf12dtheta2 = (dSf12dPf12*dPf12dtheta2+dSf12dQf12*dQf12dtheta2)*multiplier;

dSt12dVm1 =    (dSt12dPt12*dPt12dVm1+   dSt12dQt12*dQt12dVm1)*multiplier;
dSt12dVm2 =    (dSt12dPt12*dPt12dVm2+   dSt12dQt12*dQt12dVm2)*multiplier;
dSt12dtheta1 = (dSt12dPt12*dPt12dtheta1+dSt12dQt12*dQt12dtheta1)*multiplier;
dSt12dtheta2 = (dSt12dPt12*dPt12dtheta2+dSt12dQt12*dQt12dtheta2)*multiplier;


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

Sf23 = Pf23*Pf23 + Qf23*Qf23;
St23 = Pt23*Pt23 + Qt23*Qt23;

dSf23dPf23 = 2*Pf23;
dSf23dQf23 = 2*Qf23;
dSt23dPt23 = 2*Pt23;
dSt23dQt23 = 2*Qt23;

dSf23dVm2 =    (dSf23dPf23*dPf23dVm2+   dSf23dQf23*dQf23dVm2)*multiplier;
dSf23dVm3 =    (dSf23dPf23*dPf23dVm3+   dSf23dQf23*dQf23dVm3)*multiplier;
dSf23dtheta2 = (dSf23dPf23*dPf23dtheta2+dSf23dQf23*dQf23dtheta2)*multiplier;
dSf23dtheta3 = (dSf23dPf23*dPf23dtheta3+dSf23dQf23*dQf23dtheta3)*multiplier;

dSt23dVm2 =    (dSt23dPt23*dPt23dVm2+   dSt23dQt23*dQt23dVm2)*multiplier;
dSt23dVm3 =    (dSt23dPt23*dPt23dVm3+   dSt23dQt23*dQt23dVm3)*multiplier;
dSt23dtheta2 = (dSt23dPt23*dPt23dtheta2+dSt23dQt23*dQt23dtheta2)*multiplier;
dSt23dtheta3 = (dSt23dPt23*dPt23dtheta3+dSt23dQt23*dQt23dtheta3)*multiplier;

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
dPt24dtheta4 = Vm2*Vm4*(-Gtf*sin(theta42)+Btf*cos(theta42));

dQt24dVm2 =Vm4*(-Btf*cos(theta42)+Gtf*sin(theta42));
dQt24dVm4 =-2*Btt*Vm4+Vm2*(-Btf*cos(theta42)+Gtf*sin(theta42));
dQt24dtheta2 =Vm2*Vm4*(-Btf*sin(theta42)-Gtf*cos(theta42));
dQt24dtheta4 =Vm2*Vm4*(Btf*sin(theta42)+Gtf*cos(theta42));

Sf24 = Pf24*Pf24 + Qf24*Qf24;
St24 = Pt24*Pt24 + Qt24*Qt24;

dSf24dPf24 = 2*Pf24;
dSf24dQf24 = 2*Qf24;
dSt24dPt24 = 2*Pt24;
dSt24dQt24 = 2*Qt24;

dSf24dVm2 =    (dSf24dPf24*dPf24dVm2+   dSf24dQf24*dQf24dVm2)*multiplier;
dSf24dVm4 =    (dSf24dPf24*dPf24dVm4+   dSf24dQf24*dQf24dVm4)*multiplier;
dSf24dtheta2 = (dSf24dPf24*dPf24dtheta2+dSf24dQf24*dQf24dtheta2)*multiplier;
dSf24dtheta4 = (dSf24dPf24*dPf24dtheta4+dSf24dQf24*dQf24dtheta4)*multiplier;

dSt24dVm2 =    (dSt24dPt24*dPt24dVm2+   dSt24dQt24*dQt24dVm2)*multiplier;
dSt24dVm4 =    (dSt24dPt24*dPt24dVm4+   dSt24dQt24*dQt24dVm4)*multiplier;
dSt24dtheta2 = (dSt24dPt24*dPt24dtheta2+dSt24dQt24*dQt24dtheta2)*multiplier;
dSt24dtheta4 = (dSt24dPt24*dPt24dtheta4+dSt24dQt24*dQt24dtheta4)*multiplier;

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

Sf45 = Pf45*Pf45 + Qf45*Qf45;
St45 = Pt45*Pt45 + Qt45*Qt45;

dSf45dPf45 = 4*Pf45;
dSf45dQf45 = 4*Qf45;
dSt45dPt45 = 4*Pt45;
dSt45dQt45 = 4*Qt45;

dSf45dVm4 =    (dSf45dPf45*dPf45dVm4+   dSf45dQf45*dQf45dVm4)*multiplier;
dSf45dVm5 =    (dSf45dPf45*dPf45dVm5+   dSf45dQf45*dQf45dVm5)*multiplier;
dSf45dtheta4 = (dSf45dPf45*dPf45dtheta4+dSf45dQf45*dQf45dtheta4)*multiplier;
dSf45dtheta5 = (dSf45dPf45*dPf45dtheta5+dSf45dQf45*dQf45dtheta5)*multiplier;

dSt45dVm4 =    (dSt45dPt45*dPt45dVm4+   dSt45dQt45*dQt45dVm4)*multiplier;
dSt45dVm5 =    (dSt45dPt45*dPt45dVm5+   dSt45dQt45*dQt45dVm5)*multiplier;
dSt45dtheta4 = (dSt45dPt45*dPt45dtheta4+dSt45dQt45*dQt45dtheta4)*multiplier;
dSt45dtheta5 = (dSt45dPt45*dPt45dtheta5+dSt45dQt45*dQt45dtheta5)*multiplier;

% number of buses
N = 5;
% number of generators
Gen = 1;
% number of lines
Line = 4;

J=zeros(10,12);

J(3,1) = dSf12dtheta1;
J(4,1) = dSt12dtheta1;

J(3,2) = dSf12dVm1;
J(4,2) = dSt12dVm1;

J(3,3) = dSf12dtheta2;
J(4,3) = dSt12dtheta2;
J(5,3) = dSf23dtheta2;
J(6,3) = dSt23dtheta2;
J(7,3) = dSf24dtheta2;
J(8,3) = dSt24dtheta2;

J(3,4) = dSf12dVm2;
J(4,4) = dSt12dVm2;
J(5,4) = dSf23dVm2;
J(6,4) = dSt23dVm2;
J(7,4) = dSf24dVm2;
J(8,4) = dSt24dVm2;

J(5,5) = dSf23dtheta3;
J(6,5) = dSt23dtheta3;

J(1,6) = dIEC1dVm3;
J(2,6) = dIEC2dVm3;
J(5,6) = dSf23dVm3;
J(6,6) = dSt23dVm3;

J(1,8) = dIEC1dQg;
J(2,8) = dIEC2dQg;

J(7,9) = dSf24dtheta4;
J(8,9) = dSt24dtheta4;
J(9,9) = dSf45dtheta4;
J(10,9) = dSt45dtheta4;

J(7,10) = dSf24dVm4;
J(8,10) = dSt24dVm4;
J(9,10) = dSf45dVm4;
J(10,10) = dSt45dVm4;

J(9,11) = dSf45dtheta5;
J(10,11) = dSt45dtheta5;

J(9,12) = dSf45dVm5;
J(10,12) = dSt45dVm5;

J

JW(1,1) =  2*Line+2*Gen;
JW(1,2) = 2*N+2*Gen;
count = 1;

for i=1:(2*Line+2*Gen)
    for j=1:(2*N+2*Gen)
        if J(i,j)~=0
            count = count+1;
            JW(count, 1) = i;
            JW(count, 2) = j;
            JW(count, 3) = J(i,j);
        end
    end
end
writematrix(JW,'cicj.csv') 
            






