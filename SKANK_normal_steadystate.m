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
B=rB*(T);
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

%% end own model equations

for iter = 1:length(M_.params) %update parameters set in the file
  eval([ 'M_.params(' num2str(iter) ') = ' M_.param_names(iter,:) ';' ])
end

NumberOfEndogenousVariables = M_.orig_endo_nbr; %auxiliary variables are set automatically
for ii = 1:NumberOfEndogenousVariables
  varname = deblank(M_.endo_names(ii,:));
  eval(['ys(' int2str(ii) ') = ' varname ';']);
end
