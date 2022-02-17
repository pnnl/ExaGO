Vm=2;
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

Pf12=Gff*Vm*Vm+Vm*Vm*(Gft*cos(theta12)+Bft*sin(theta12));
Qf12=-Bff*Vm*Vm+Vm*Vm*(-Bft*cos(theta12)+Gft*sin(theta12));
Pt12=Gtt*Vm*Vm+Vm*Vm*(Gtf*cos(theta21)+Btf*sin(theta21));
Qt12=-Btt*Vm*Vm+Vm*Vm*(-Btf*cos(theta21)+Gtf*sin(theta21));

Sf12 = Pf12*Pf12 + Qf12*Qf12;
St12 = Pt12*Pt12 + Qt12*Qt12;

Pf23=Gfft*Vm*Vm+Vm*Vm*(Gftt*cos(theta23)+Bftt*sin(theta23));
Qf23=-Bfft*Vm*Vm+Vm*Vm*(-Bftt*cos(theta23)+Gftt*sin(theta23));
Pt23=Gttt*Vm*Vm+Vm*Vm*(Gtft*cos(theta32)+Btft*sin(theta32));
Qt23=-Bttt*Vm*Vm+Vm*Vm*(-Btft*cos(theta32)+Gtft*sin(theta32));

Sf23 = Pf23*Pf23 + Qf23*Qf23;
St23 = Pt23*Pt23 + Qt23*Qt23;

Pf24=Gff*Vm*Vm+Vm*Vm*(Gft*cos(theta24)+Bft*sin(theta24));
Qf24=-Bff*Vm*Vm+Vm*Vm*(-Bft*cos(theta24)+Gft*sin(theta24));
Pt24=Gtt*Vm*Vm+Vm*Vm*(Gtf*cos(theta42)+Btf*sin(theta42));
Qt24=-Btt*Vm*Vm+Vm*Vm*(-Btf*cos(theta42)+Gtf*sin(theta42));

Sf24 = Pf24*Pf24 + Qf24*Qf24;
St24 = Pt24*Pt24 + Qt24*Qt24;

Pf45=Gff*Vm*Vm+Vm*Vm*(Gft*cos(theta45)+Bft*sin(theta45));
Qf45=-Bff*Vm*Vm+Vm*Vm*(-Bft*cos(theta45)+Gft*sin(theta45));
Pt45=Gtt*Vm*Vm+Vm*Vm*(Gtf*cos(theta54)+Btf*sin(theta54));
Qt45=-Btt*Vm*Vm+Vm*Vm*(-Btf*cos(theta54)+Gtf*sin(theta54));

Sf45 = Pf45*Pf45 + Qf45*Qf45;
St45 = Pt45*Pt45 + Qt45*Qt45;


Res(1)=IEC1;
Res(2)=IEC2;
Res(3)=Sf12;
Res(4)=St12;
Res(5)=Sf23;
Res(6)=St23;
Res(7)=Sf24;
Res(8)=St24;
Res(9)=Sf45;
Res(10)=St45;

Res

ResW(1,1) =  10;

for i=1:10
    ResW(i+1,1) = Res(i);
end
writematrix(ResW,'cic.csv') 


