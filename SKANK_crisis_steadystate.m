function [ys,check] = SKANK_normal_steadystate(ys,exe)
global M_

% read out parameters to access them with their name
NumberOfParameters = M_.param_nbr;
for ii = 1:NumberOfParameters
  paramname = deblank(M_.param_names(ii,:));
  eval([ paramname ' = M_.params(' int2str(ii) ');']);
end
% initialize indicator
check = 0;

%% Enter model equations here

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

siggma=siggma_ss;
piee=piee_ss;
bR=bR_ss; 
D_ss=D; 
R_ss=R; 
Y_ss=Y; 
F_ss=F;

%% end own model equations

for iter = 1:length(M_.params) %update parameters set in the file
  eval([ 'M_.params(' num2str(iter) ') = ' M_.param_names(iter,:) ';' ])
end

NumberOfEndogenousVariables = M_.orig_endo_nbr; %auxiliary variables are set automatically
for ii = 1:NumberOfEndogenousVariables
  varname = deblank(M_.endo_names(ii,:));
  eval(['ys(' int2str(ii) ') = ' varname ';']);
end
