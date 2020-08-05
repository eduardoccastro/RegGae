function [SS,EP] = RegGae_steadystate(Gamma)
%% RegGae: a toolkit for macroprudencial policy
% Copyright 2020, Eduardo C. Castro ecastro@bcb.gov.br
% May be used at your own risk but proper credit is required
% Function computes the SS based on the params vector Gamma for the crisis model slack 

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
bR_ss=Gamma(38);

%Enter the process that finds the SS:
RD=piee_ss/betta_ss/(1-zetaD);
RF=piee_ss/betta_ss/(1-zetaF);
RN=piee_ss/betta_ss/(1-zetaD);
R=piee_ss/betta_ss;
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

%Equations with double loss and with lump-sum tax
T =(betta_ss*(F - M)*(RN*piee_ss - RD*piee_ss - 2*RD*RN*zetaD + RD*piee_ss*zetaD + RD*RN*zetaD^2))/(RD*piee_ss^2*rB - RN*piee_ss^2*rB + RD*RN*betta_ss*zetaD^2 - RD*piee_ss^2*rB*zetaD - RD*RT*betta_ss*rB + RN*RT*betta_ss*rB - 2*RD*RN*betta_ss*zetaD + RD*RT*betta_ss*zetaD + RN*RT*betta_ss*zetaD + 2*RD*RN*betta_ss*rB*zetaD - RN*RT*betta_ss*rB*zetaD - RD*RN*betta_ss*rB*zetaD^2) - (betta_ss*(M + F*RF*(zetaF - 1))*(RD - RN + RN*zetaD))/(RD*piee_ss^2*rB - RN*piee_ss^2*rB + RD*RN*betta_ss*zetaD^2 - RD*piee_ss^2*rB*zetaD - RD*RT*betta_ss*rB + RN*RT*betta_ss*rB - 2*RD*RN*betta_ss*zetaD + RD*RT*betta_ss*zetaD + RN*RT*betta_ss*zetaD + 2*RD*RN*betta_ss*rB*zetaD - RN*RT*betta_ss*rB*zetaD - RD*RN*betta_ss*rB*zetaD^2) - (betta_ss*piee_ss*(OmegaF + M*(1/piee_ss - 1))*(RN - RD + RD*zetaD))/(RD*piee_ss^2*rB - RN*piee_ss^2*rB + RD*RN*betta_ss*zetaD^2 - RD*piee_ss^2*rB*zetaD - RD*RT*betta_ss*rB + RN*RT*betta_ss*rB - 2*RD*RN*betta_ss*zetaD + RD*RT*betta_ss*zetaD + RN*RT*betta_ss*zetaD + 2*RD*RN*betta_ss*rB*zetaD - RN*RT*betta_ss*rB*zetaD - RD*RN*betta_ss*rB*zetaD^2);
D =((F - M)*(RN*betta_ss*piee_ss - RT*betta_ss*piee_ss + RN*piee_ss^2*rB - RN*RT*betta_ss*rB - RN*RT*betta_ss*zetaD - RN*betta_ss*piee_ss*rB + RT*betta_ss*piee_ss*rB + RN*RT*betta_ss*rB*zetaD))/(RD*piee_ss^2*rB - RN*piee_ss^2*rB + RD*RN*betta_ss*zetaD^2 - RD*piee_ss^2*rB*zetaD - RD*RT*betta_ss*rB + RN*RT*betta_ss*rB - 2*RD*RN*betta_ss*zetaD + RD*RT*betta_ss*zetaD + RN*RT*betta_ss*zetaD + 2*RD*RN*betta_ss*rB*zetaD - RN*RT*betta_ss*rB*zetaD - RD*RN*betta_ss*rB*zetaD^2) + ((M + F*RF*(zetaF - 1))*(rB*piee_ss^2 + RN*betta_ss - RT*betta_ss - RN*betta_ss*rB - RN*betta_ss*zetaD + RN*betta_ss*rB*zetaD))/(RD*piee_ss^2*rB - RN*piee_ss^2*rB + RD*RN*betta_ss*zetaD^2 - RD*piee_ss^2*rB*zetaD - RD*RT*betta_ss*rB + RN*RT*betta_ss*rB - 2*RD*RN*betta_ss*zetaD + RD*RT*betta_ss*zetaD + RN*RT*betta_ss*zetaD + 2*RD*RN*betta_ss*rB*zetaD - RN*RT*betta_ss*rB*zetaD - RD*RN*betta_ss*rB*zetaD^2) - (betta_ss*piee_ss*(OmegaF + M*(1/piee_ss - 1))*(RN - RT - RN*rB + RT*rB))/(RD*piee_ss^2*rB - RN*piee_ss^2*rB + RD*RN*betta_ss*zetaD^2 - RD*piee_ss^2*rB*zetaD - RD*RT*betta_ss*rB + RN*RT*betta_ss*rB - 2*RD*RN*betta_ss*zetaD + RD*RT*betta_ss*zetaD + RN*RT*betta_ss*zetaD + 2*RD*RN*betta_ss*rB*zetaD - RN*RT*betta_ss*rB*zetaD - RD*RN*betta_ss*rB*zetaD^2);
N =((M + F*RF*(zetaF - 1))*(rB*piee_ss^2 + RD*betta_ss - RT*betta_ss - RD*betta_ss*rB))/(RD*piee_ss^2*rB - RN*piee_ss^2*rB + RD*RN*betta_ss*zetaD^2 - RD*piee_ss^2*rB*zetaD - RD*RT*betta_ss*rB + RN*RT*betta_ss*rB - 2*RD*RN*betta_ss*zetaD + RD*RT*betta_ss*zetaD + RN*RT*betta_ss*zetaD + 2*RD*RN*betta_ss*rB*zetaD - RN*RT*betta_ss*rB*zetaD - RD*RN*betta_ss*rB*zetaD^2) + ((F - M)*(RD*betta_ss*piee_ss - RT*betta_ss*piee_ss + RD*piee_ss^2*rB - RD*piee_ss^2*rB*zetaD - RD*RT*betta_ss*rB + RD*RT*betta_ss*zetaD - RD*betta_ss*piee_ss*rB + RT*betta_ss*piee_ss*rB - RD*betta_ss*piee_ss*zetaD + RD*betta_ss*piee_ss*rB*zetaD))/(RD*piee_ss^2*rB - RN*piee_ss^2*rB + RD*RN*betta_ss*zetaD^2 - RD*piee_ss^2*rB*zetaD - RD*RT*betta_ss*rB + RN*RT*betta_ss*rB - 2*RD*RN*betta_ss*zetaD + RD*RT*betta_ss*zetaD + RN*RT*betta_ss*zetaD + 2*RD*RN*betta_ss*rB*zetaD - RN*RT*betta_ss*rB*zetaD - RD*RN*betta_ss*rB*zetaD^2) - (betta_ss*piee_ss*(OmegaF + M*(1/piee_ss - 1))*(RD - RT - RD*rB + RT*rB - RD*zetaD + RD*rB*zetaD))/(RD*piee_ss^2*rB - RN*piee_ss^2*rB + RD*RN*betta_ss*zetaD^2 - RD*piee_ss^2*rB*zetaD - RD*RT*betta_ss*rB + RN*RT*betta_ss*rB - 2*RD*RN*betta_ss*zetaD + RD*RT*betta_ss*zetaD + RN*RT*betta_ss*zetaD + 2*RD*RN*betta_ss*rB*zetaD - RN*RT*betta_ss*rB*zetaD - RD*RN*betta_ss*rB*zetaD^2);

B=rB*T; 
b=N/(wD*D*(1-zetaD)+wF*F*(1-zetaF));
bI=RN*N/(wD*(1-zetaD)*RD*D+wF*(1-zetaF)*RF*F);
OmegaD=zetaD*RD*D/piee_ss;
betta=betta_ss;
A=A_ss;
bR=bR_ss; 
siggma=siggma_ss;
piee=piee_ss;

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

EP(1,1)=bR_ss;
EP(2,1)=D; 
EP(3,1)=R; 
EP(4,1)=F;
EP(5,1)=Y; 


end

