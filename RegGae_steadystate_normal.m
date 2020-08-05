function [SS,EP] = RegGae_steadystate(Gamma)
%% RegGae: a toolkit for macroprudencial policy
% Copyright 2020, Eduardo C. Castro ecastro@bcb.gov.br
% May be used at your own risk but proper credit is required
% Function computes the SS based on the params vector Gamma for the normal model tight.

%Enter constant parameters here:
s=Gamma(1); 
a=Gamma(2); 
zetaF=Gamma(3);
zetaD=Gamma(4);
wF=Gamma(5); 
wD=Gamma(6);
m=Gamma(7); 
delta=Gamma(8); 
phi=Gamma(9); 
theta=Gamma(10); 
alfa=Gamma(11); 
siggma_ss=Gamma(12); 
betta_ss=Gamma(13);
rho_A=Gamma(14); 
rho_R=Gamma(15); 
rho_bR=Gamma(16); 
rho_beta=Gamma(17); 
rho_sigma=Gamma(18); 
siggma_sigma=Gamma(19);
siggma_A=Gamma(20); 
siggma_R=Gamma(21); 
siggma_bR=Gamma(22); 
siggma_pi=Gamma(23); 
siggma_beta=Gamma(24);
epps=Gamma(25); 
phi_pi=Gamma(26); 
phi_y=Gamma(27); 
phi_DY=Gamma(28); 
csi=Gamma(29); 
vB=Gamma(30); 
A_ss=Gamma(31); 
piee_ss=Gamma(32); 
b_ss=Gamma(33); 
nu=Gamma(34); 
Q=Gamma(35); 
rB=Gamma(36); 
h=Gamma(37);

%Enter the process that finds the SS:
RD=piee_ss/betta_ss/(1-zetaD);
RF=piee_ss/betta_ss/(1-zetaF);
RN=piee_ss/betta_ss/(1-zetaD);
RT=piee_ss/betta_ss;
 
rK=(1-nu)/betta_ss+(1-zetaF)*RF/piee_ss*nu-(1-delta-zetaF);
piee_i=((piee_ss^(1-epps)-csi)/(1-csi))^(1/(1-epps));
Pd=(1-csi)/(1-csi*piee_ss^epps)*(piee_ss/piee_i)^epps;
SRi=(piee_i/piee_ss)*((epps-1)/epps);
SR=SRi*Pd;
Y_K=rK/(alfa*SR);
K_L=(Y_K*Pd/A_ss)^(1/(alfa-1));
I_L=(delta+zetaF)*K_L;
C_L=Y_K*K_L-I_L;
W=(((SR*A_ss)/Pd)*(alfa/rK)^alfa)^(1/(1-alfa))*(1-alfa);

L=(((1-h)*C_L)^(-siggma_ss)*W/(theta*(1+m-betta_ss*m/piee_ss)))^(1/(phi+siggma_ss));
K=K_L*L;
Y=A_ss/Pd*K^alfa*L^(1-alfa);
I=(delta+zetaF)*K;
X1=SRi*Y*piee_ss^(1+epps)/(1-csi*betta_ss);
X2=piee_ss^epps*Y/(1-csi*betta_ss);
C=C_L*L;
Lambda=theta/W*L^phi;
OmegaS=C-W*L;
M=m*C;
F=nu*K;
OmegaF=rK*K-I+F-(1-zetaF)*RF/piee_ss*F;

R=1.01;  
erro=1;

while 0.00001<erro | -0.00001>erro

T =(piee_ss*(OmegaF + M*(1/piee_ss - 1))*(b_ss*wD*zetaD - b_ss*wD + 1))/(RD - RT - RD*rB - RN*b_ss*wD + RT*b_ss*wD + R*piee_ss*rB + 2*RN*b_ss*wD*zetaD - RT*b_ss*wD*zetaD - RN*b_ss*wD*zetaD^2 + RN*b_ss*rB*wD - R*b_ss*piee_ss*rB*wD - 2*RN*b_ss*rB*wD*zetaD + RN*b_ss*rB*wD*zetaD^2 + R*b_ss*piee_ss*rB*wD*zetaD) - ((F - M)*(piee_ss - RD + RN*b_ss*wD - b_ss*piee_ss*wD - 2*RN*b_ss*wD*zetaD + b_ss*piee_ss*wD*zetaD + RN*b_ss*wD*zetaD^2))/(RD - RT - RD*rB - RN*b_ss*wD + RT*b_ss*wD + R*piee_ss*rB + 2*RN*b_ss*wD*zetaD - RT*b_ss*wD*zetaD - RN*b_ss*wD*zetaD^2 + RN*b_ss*rB*wD - R*b_ss*piee_ss*rB*wD - 2*RN*b_ss*rB*wD*zetaD + RN*b_ss*rB*wD*zetaD^2 + R*b_ss*piee_ss*rB*wD*zetaD) + (F*b_ss*wF*(zetaF - 1)*(RD - RN + RN*zetaD))/(RD - RT - RD*rB - RN*b_ss*wD + RT*b_ss*wD + R*piee_ss*rB + 2*RN*b_ss*wD*zetaD - RT*b_ss*wD*zetaD - RN*b_ss*wD*zetaD^2 + RN*b_ss*rB*wD - R*b_ss*piee_ss*rB*wD - 2*RN*b_ss*rB*wD*zetaD + RN*b_ss*rB*wD*zetaD^2 + R*b_ss*piee_ss*rB*wD*zetaD); 
D =((F - M)*(RT - piee_ss + piee_ss*rB - R*piee_ss*rB))/(RD - RT - RD*rB - RN*b_ss*wD + RT*b_ss*wD + R*piee_ss*rB + 2*RN*b_ss*wD*zetaD - RT*b_ss*wD*zetaD - RN*b_ss*wD*zetaD^2 + RN*b_ss*rB*wD - R*b_ss*piee_ss*rB*wD - 2*RN*b_ss*rB*wD*zetaD + RN*b_ss*rB*wD*zetaD^2 + R*b_ss*piee_ss*rB*wD*zetaD) - (piee_ss*(OmegaF + M*(1/piee_ss - 1))*(rB - 1))/(RD - RT - RD*rB - RN*b_ss*wD + RT*b_ss*wD + R*piee_ss*rB + 2*RN*b_ss*wD*zetaD - RT*b_ss*wD*zetaD - RN*b_ss*wD*zetaD^2 + RN*b_ss*rB*wD - R*b_ss*piee_ss*rB*wD - 2*RN*b_ss*rB*wD*zetaD + RN*b_ss*rB*wD*zetaD^2 + R*b_ss*piee_ss*rB*wD*zetaD) - (F*b_ss*wF*(zetaF - 1)*(RN - RT - RN*rB - RN*zetaD + R*piee_ss*rB + RN*rB*zetaD))/(RD - RT - RD*rB - RN*b_ss*wD + RT*b_ss*wD + R*piee_ss*rB + 2*RN*b_ss*wD*zetaD - RT*b_ss*wD*zetaD - RN*b_ss*wD*zetaD^2 + RN*b_ss*rB*wD - R*b_ss*piee_ss*rB*wD - 2*RN*b_ss*rB*wD*zetaD + RN*b_ss*rB*wD*zetaD^2 + R*b_ss*piee_ss*rB*wD*zetaD);
N =(b_ss*piee_ss*wD*(OmegaF + M*(1/piee_ss - 1))*(rB - 1)*(zetaD - 1))/(RD - RT - RD*rB - RN*b_ss*wD + RT*b_ss*wD + R*piee_ss*rB + 2*RN*b_ss*wD*zetaD - RT*b_ss*wD*zetaD - RN*b_ss*wD*zetaD^2 + RN*b_ss*rB*wD - R*b_ss*piee_ss*rB*wD - 2*RN*b_ss*rB*wD*zetaD + RN*b_ss*rB*wD*zetaD^2 + R*b_ss*piee_ss*rB*wD*zetaD) - (b_ss*wD*(F - M)*(zetaD - 1)*(RT - piee_ss + piee_ss*rB - R*piee_ss*rB))/(RD - RT - RD*rB - RN*b_ss*wD + RT*b_ss*wD + R*piee_ss*rB + 2*RN*b_ss*wD*zetaD - RT*b_ss*wD*zetaD - RN*b_ss*wD*zetaD^2 + RN*b_ss*rB*wD - R*b_ss*piee_ss*rB*wD - 2*RN*b_ss*rB*wD*zetaD + RN*b_ss*rB*wD*zetaD^2 + R*b_ss*piee_ss*rB*wD*zetaD) - (F*b_ss*wF*(zetaF - 1)*(RD - RT - RD*rB + R*piee_ss*rB))/(RD - RT - RD*rB - RN*b_ss*wD + RT*b_ss*wD + R*piee_ss*rB + 2*RN*b_ss*wD*zetaD - RT*b_ss*wD*zetaD - RN*b_ss*wD*zetaD^2 + RN*b_ss*rB*wD - R*b_ss*piee_ss*rB*wD - 2*RN*b_ss*rB*wD*zetaD + RN*b_ss*rB*wD*zetaD^2 + R*b_ss*piee_ss*rB*wD*zetaD);
B=rB*T;
bI=RN*N/(wD*(1-zetaD)*RD*D+wF*(1-zetaF)*RF*F);
OmegaD=zetaD*RD/piee_ss*D;
bR=2*b_ss-vB-bI;
bR_ss=bR;
b=b_ss;
betta=betta_ss;
A=A_ss;
siggma=siggma_ss;
piee=piee_ss;
D_ss=D;
F_ss=F;
R_ss=R;
Y_ss=Y;
erro=R*B+(1-zetaD)*RD*D+(1-zetaF)*RF*F-m*C-RT*T-RN*N;
R=R-erro/10;

end

%Define the SS output (in the same order of the .mod file)
SS=[];
SS(1,1)=Y;
SS(2,1)=I;
SS(3,1)=L;
SS(4,1)=W;
SS(5,1)=b;
SS(6,1)=bI;
SS(7,1)=piee_i;
SS(8,1)=SR;
SS(9,1)=SRi;
SS(10,1)=OmegaS;
SS(11,1)=OmegaF;
SS(12,1)=OmegaD;
SS(13,1)=B;
SS(14,1)=T;
SS(15,1)=F;
SS(16,1)=D;
SS(17,1)=N;
SS(18,1)=R;
SS(19,1)=RF;
SS(20,1)=RD;
SS(21,1)=RT;
SS(22,1)=betta;
SS(23,1)=C;
SS(24,1)=Pd;
SS(25,1)=A;
SS(26,1)=bR;
SS(27,1)=K;
SS(28,1)=siggma;
SS(29,1)=X1;
SS(30,1)=X2;
SS(31,1)=piee;
SS(32,1)=Lambda;
SS(33,1)=RN;
SS(34,1)=rK;

% EP: endogenous parameters
EP(1,1)=bR_ss;
EP(2,1)=D_ss; 
EP(3,1)=R_ss; 
EP(4,1)=F_ss;
EP(5,1)=Y_ss; 

end

