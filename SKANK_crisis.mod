// SKANK: Structural Capital Adequacy in a NK Model
// Copyright 2020 Eduardo Castro. All rights reserved.
// Central Bank of Brazil ecastro@bcb.gov.br

// Crisis model: macroprudencial rule is slack

var  // Total: 
// Ordering must be respected: 
// First: static variables 
    Y I L W b bI piee_i SR SRi OmegaS OmegaF OmegaD 

// Second: Pure predetermined variables
  B T F D N R RF RD RT betta C Pd A bR K siggma 

// Third: mixed (lead-lag) variables

// Fourth: Pure forward-looking variables 
    X1 X2 piee Lambda RN rK;

// Alert: variables in all regime-.mod files must be declared in the same order

varexo eps_A eps_R eps_bR eps_pi eps_beta eps_sigma;

parameters //exogenous
s a zetaF zetaD wF wD m delta phi theta alfa siggma_ss betta_ss
rho_A rho_R rho_bR rho_beta rho_sigma siggma_sigma
siggma_A siggma_R siggma_bR siggma_pi siggma_beta
epps phi_pi phi_y phi_DY csi vB A_ss piee_ss b_ss nu Q rB h

//endogenous 
bR_ss D_ss R_ss F_ss Y_ss ; //Variable parameter should be declared last

set_param_value('s',Gamma(1,Reg));
set_param_value('a',Gamma(2,Reg));
set_param_value('zetaF',Gamma(3,Reg));
set_param_value('zetaD',Gamma(4,Reg)); 
set_param_value('wF',Gamma(5,Reg)); 
set_param_value('wD',Gamma(6,Reg)); 
set_param_value('m',Gamma(7,Reg)); 
set_param_value('delta',Gamma(8,Reg));
set_param_value('phi',Gamma(9,Reg)); 
set_param_value('theta',Gamma(10,Reg));
set_param_value('alfa',Gamma(11,Reg)); 
set_param_value('siggma_ss',Gamma(12,Reg));
set_param_value('betta_ss',Gamma(13,Reg));
set_param_value('rho_A',Gamma(14,Reg));
set_param_value('rho_R',Gamma(15,Reg)); 
set_param_value('rho_bR',Gamma(16,Reg));
set_param_value('rho_beta',Gamma(17,Reg));
set_param_value('rho_sigma',Gamma(18,Reg));
set_param_value('siggma_sigma',Gamma(19,Reg));
set_param_value('siggma_A',Gamma(20,Reg));
set_param_value('siggma_R',Gamma(21,Reg));
set_param_value('siggma_bR',Gamma(22,Reg));
set_param_value('siggma_pi',Gamma(23,Reg));
set_param_value('siggma_beta',Gamma(24,Reg));
set_param_value('epps',Gamma(25,Reg));
set_param_value('phi_pi',Gamma(26,Reg));
set_param_value('phi_y',Gamma(27,Reg));
set_param_value('phi_DY',Gamma(28,Reg));
set_param_value('csi',Gamma(29,Reg));
set_param_value('vB',Gamma(30,Reg));
set_param_value('A_ss',Gamma(31,Reg));
set_param_value('piee_ss',Gamma(32,Reg));
set_param_value('b_ss',Gamma(33,Reg)); 
set_param_value('nu',Gamma(34,Reg));
set_param_value('Q',Gamma(35,Reg));
set_param_value('rB',Gamma(36,Reg));
set_param_value('h',Gamma(37,Reg));
set_param_value('bR_ss',Gamma(38,Reg));

model;
B+D+F=m*C+T+N; 
RN*N(-1)=R(-1)*B(-1)+(1-zetaD)*RD(-1)*D(-1)+(1-zetaF)*RF(-1)*F(-1)-m*C(-1)-RT(-1)*T(-1); 
N=b*(wD*D*(1-zetaD)+wF*F*(1-zetaF)); 
R=RD*(1-zetaD);  
R=RF*(1-zetaF);
R=RT;
B=rB*T;
bI=RN*N(-1)/(wD*(1-zetaD)*RD(-1)*D(-1)+wF*(1-zetaF)*RF(-1)*F(-1)); 
//bR=bR_ss+siggma_bR*eps_bR;
bR=bR(-1)*rho_bR+(1-rho_bR)*bR_ss*((D+F)/(D_ss+F_ss))^phi_DY+siggma_bR*eps_bR; 
//1+bR=(1+bR_ss)*(bR(-1)/bR_ss)^rho_bR*(((D+F)/(D_ss+F_ss))^phi_DY)^(1-rho_bR)*exp(siggma_bR*eps_bR);
OmegaD=zetaD*D(-1)*RD(-1)/piee;
(1+m)*C+(1-zetaD)*RD(-1)*D(-1)/piee+N+T+OmegaD+(R-1)*B=W*L+D+m*C(-1)/piee+OmegaS+OmegaF+(1-zetaD)*RN*N(-1)/piee+RT(-1)/piee*T(-1);
(C-h*C(-1))^(-siggma)-Lambda*(1+m)+betta*Lambda(+1)*m/piee(+1)=0;
theta*L^(phi)=Lambda*W;
Lambda=betta*(1-zetaD)*Lambda(+1)*RN(+1)/piee(+1); 
Lambda=betta*(1-zetaD)*Lambda(+1)*RD/piee(+1);
//Lambda=betta*Lambda(+1)*RT/piee(+1); 
Y=(A/Pd)*K(-1)^alfa*L^(1-alfa);
ln(A)=rho_A*ln(A(-1))+(1-rho_A)*ln(A_ss)+siggma_A*eps_A;
W=(1-alfa)*SR*Y/L;
rK=alfa*SR*Y/K(-1);
SR=SRi*Pd;
piee_i=(epps/(epps-1))*(X1/X2);
X1=SRi*piee^(1+epps)*Y+csi*betta*(Lambda(+1)/Lambda)*X1(+1);
X2=(piee^(epps))*Y+csi*betta*Lambda(+1)/Lambda*X2(+1);
piee=(((1-csi)*piee_i^(1-epps)+csi)^(1/(1-epps)))+siggma_pi*eps_pi;
Pd=(1-csi)*(piee/piee_i)^epps+csi*piee^epps*Pd(-1);
OmegaF=(1+rK-delta-zetaF)*K(-1)-K+F-(1-zetaF)*F(-1)*RF(-1)/piee;
K=I+(1-delta-zetaF)*K(-1); 
nu*K=F;
Lambda*(nu-Q)+betta*Lambda(+1)*(rK(+1)- (1-zetaF)*RF/piee(+1)*nu+Q(+1)*(1-delta-zetaF))=0;
R=R_ss*(R(-1)/R_ss)^(rho_R)*((piee/piee_ss)^phi_pi*(Y/Y_ss)^phi_y)^(1-rho_R)*exp(siggma_R*eps_R);
OmegaS=C-W*L; 
Y=C+I;
ln(siggma)=rho_sigma*ln(siggma(-1))+(1-rho_sigma)*ln(siggma_ss)+siggma_sigma*eps_sigma;
betta=rho_beta*betta(-1)+(1-rho_beta)*betta_ss-siggma_beta*eps_beta;
end;

shocks;
var eps_A=1;
//var eps_R=1;
//var eps_bR=1;
//var eps_beta=1;
//var eps_pi=1;
//var eps_sigma=1;
end;

steady;
check;
//stoch_simul(order=1, noprint, irf=30); // , 
// Options irf_shocks=(eps_A, eps_R, eps_bR, eps_beta, eps_pi) 