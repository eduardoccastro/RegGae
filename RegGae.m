%% RegGae: a toolkit for macroprudencial policy
% Copyright 2020, Eduardo C. Castro ecastro@bcb.gov.br
% May be used at your own risk but proper credit is required

% This file is setup to implement RegGae for two dynare files for the DSGE
% SKANK (forthcoming article) with two regimes, normal (tight) and crisis
% (slack).

clc; addpath c:\dynare\4.5.7\matlab; clear 

% I - Preparation: Regime matrices (R), Transition matrices (TM) by history type iota
R=struct(); TM=struct(); % initiate struct of transition matrices
T=21; bR_ss=0.105; % in the slack model
vareps=[0;0;0;0;0;0]; % set shock for IRFs: eps_A eps_R eps_bR eps_pi eps_beta eps_sigma
Iota=[115 114 115 115 115 115 115 115 115 115 115 115 115 115 115 115 115 115 115  115 115 ];

% II - Creates matrix Gamma of parameter values of regimes
%Regimes 1-99 normal (tight) model 
Gamma=[];
Gamma(1,1)= 1;%s
Gamma(2,1)= 0.5;%a 
Gamma(3,1)= 0.015;%zetaF multiplicatin by s below
Gamma(4,1)= 0.03;%zetaD multiplicatin by s below
Gamma(5,1)= 0.5;%wF 
Gamma(6,1)= 0.8;%wD 
Gamma(7,1)= .8;%m 
Gamma(8,1)= 0.02;%delta 
Gamma(9,1)= 1;%phi 
Gamma(10,1)= 14.71;%theta 
Gamma(11,1)= 0.37;%alfa 
Gamma(12,1)= 1;%siggma_ss 
Gamma(13,1)= 0.99; %betta_ss
Gamma(14,1)= 0.75;%rho_A 
Gamma(15,1)= 0.75;%rho_R 
Gamma(16,1)= 0.25;%rho_bR 
Gamma(17,1)= 0.9;%rho_beta 
Gamma(18,1)= 0.95;%rho_sigma 
Gamma(19,1)= 0.01;%siggma_sigma
Gamma(20,1)= 0.01;%siggma_A 
Gamma(21,1)= 0.0025;%siggma_R 
Gamma(22,1)= 0.0025;%siggma_bR 
Gamma(23,1)= 0.0025;%siggma_pi 
Gamma(24,1)= 0.0025;%siggma_beta
Gamma(25,1)= 6;%epps 
Gamma(26,1)= 2.5;%phi_pi 
Gamma(27,1)= 1;%phi_y 
Gamma(28,1)= -.01;%phi_DY 
Gamma(29,1)= 0.5;%csi 
Gamma(30,1)= 0.01;%vB 
Gamma(31,1)= 1;%A_ss 
Gamma(32,1)= 1.035^(1/4);%piee_ss 
Gamma(33,1)= 0.14109;%b_ss 
Gamma(34,1)= 0.4;%nu 
Gamma(35,1)= 1;%Q 
Gamma(36,1)= 0.1;%rB 
Gamma(37,1)= 0.435;%h .435 why?
% Endogenous parameters: write 0; it will be replaced by SS values
Gamma(38,1)= 0;%bR_ss endogenous 
Gamma(39,1)= 0;%D_ss endogenous
Gamma(40,1)= 0;%R_ss endogenous
Gamma(41,1)= 0;%F_ss endogenous
Gamma(42,1)= 0;%Y_ss endogenous

for Reg=2:10
Gamma(:,Reg)=Gamma(:,1);
end
Gamma(1,1:10)=linspace(1,1,10); % socrecovery path

% Regimes number 100-onwards: crisis (slack) model
Gamma(1,101)= 1;%s
Gamma(2,101)= 0.5;%a 
Gamma(3,101)= 0.015;%zetaF multiplication by s below
Gamma(4,101)= 0.03;%zetaD multiplicatin by s below
Gamma(5,101)= 0.5;%wF 
Gamma(6,101)= 0.8;%wD 
Gamma(7,101)= 1*Gamma(7,1);%m 
Gamma(8,101)= 0.02;%delta 
Gamma(9,101)= 1;%phi 
Gamma(10,101)= 14.71;%theta 
Gamma(11,101)= 0.37;%alfa 
Gamma(12,101)= 1;%siggma_ss 
Gamma(13,101)= 0.99; %betta_ss
Gamma(14,101)= 0;%rho_A 
Gamma(15,101)= 0.75;%rho_R 
Gamma(16,101)= 0;%rho_bR 
Gamma(17,101)= 0.9;%rho_beta 
Gamma(18,101)= 0.95;%rho_sigma 
Gamma(19,101)= 0.01;%siggma_sigma
Gamma(20,101)= 0.01;%siggma_A 
Gamma(21,101)= 0.0025;%siggma_R 
Gamma(22,101)= 0.0018;%siggma_bR 
Gamma(23,101)= 0.0025;%siggma_pi 
Gamma(24,101)= 0.0025;%siggma_beta
Gamma(25,101)= 6;%epps 
Gamma(26,101)= 2.5;%phi_pi 
Gamma(27,101)= 1;%phi_y 
Gamma(28,101)= 0.0;%phi_DY 
Gamma(29,101)= 0.5;%csi 
Gamma(30,101)= 0.01;%vB 
Gamma(31,101)= 1;%A_ss 
Gamma(32,101)= 1.035^(1/4);%piee_ss 
Gamma(33,101)= 0.385780034603496;%b_ss free parameter
Gamma(34,101)= 0.4;%nu 
Gamma(35,101)= 1;%Q 
Gamma(36,101)= 0.1;%rB 
Gamma(37,101)= 0.435;%h
Gamma(38,101)= bR_ss; %%free parameter
Gamma(39,101)= 0;%D_ss endogenous
Gamma(40,101)= 0;%R_ss endogenous
Gamma(41,101)= 0;%F_ss endogenous
Gamma(42,101)= 0;%Y_ss endogenous

for Reg=101:115
Gamma(:,Reg)=Gamma(:,101);
end

Gamma(7,114)=.808;

%vardelta=.0; Gamma(1,101:115)=[1 1+vardelta/2 1+3/4*vardelta 1+7/8*vardelta linspace(1+vardelta,1,10) 1];
%{
for i=102:114
    Gamma(31,101)=.99;
Gamma(31,i)=exp(.75*log(Gamma(31,i-1)));
end
%}
%Gamma(3:4,:)=Gamma(1,:).* Gamma(3:4,:); %multiplies the NPL rates
%Gamma(33,:)=Gamma(1,:).* Gamma(33,:); %multiplies the b_ss
%Gamma(34,:)=Gamma(1,:).* Gamma(34,:); %multiplies the nu
%Gamma(10,:)= Gamma(10,:).*Gamma(1,:); %multiplies theta 
[1:25;Gamma(:,1:10) Gamma(:,101:115)];

%% III - Regime-wise linearization

%First time: execute Dynare once for each model

Reg=1; %Regime
[R.T(:,:,Reg), R.S(:,:,Reg), R.Q(:,:,Reg), R.Z(:,:,Reg), R.C(:,:,Reg), R.D(:,:,Reg), ...
    R.SSd(:,Reg), R.st(Reg), R.p(Reg), R.j(Reg), R.Ast(:,:,:,Reg), R.SSst(:,Reg), R.Dst(:,:,Reg) ]=...
    RegGae_matrix('SKANK_normal');

R.SS(:,Reg)=[R.SSst(:,Reg);R.SSd(:,Reg)];
Gamma(38:42,Reg)=M_.params(38:42); %Endogenous parameters

for i=2:10
[SS,EP]=RegGae_steadystate_normal(Gamma(:,i)); %obtains SS and EP endog. params.
Gamma(38:42,i)=EP;
[R.T(:,:,i), R.S(:,:,i), R.Q(:,:,i), R.Z(:,:,i), R.C(:,:,i), R.D(:,:,i), ...
    R.SSd(:,i), R.st(i), R.p(i), R.j(i), R.Ast(:,:,:,i), R.SSst(:,i), R.Dst(:,:,i) ]=...
    RegGae_jacobian('SKANK_normal',SS,Gamma(:,i));
R.SS(:,i)=[R.SSst(:,i); R.SSd(:,i)];
end

Reg=101;
[R.T(:,:,Reg), R.S(:,:,Reg), R.Q(:,:,Reg), R.Z(:,:,Reg), R.C(:,:,Reg), R.D(:,:,Reg), ...
    R.SSd(:,Reg), R.st(Reg), R.p(Reg), R.j(Reg), R.Ast(:,:,:,Reg), R.SSst(:,Reg), R.Dst(:,:,Reg) ]=...
    RegGae_matrix('SKANK_crisis');
Gamma(39:42,Reg)=M_.params(39:42); % endogenous parameters
R.SS(:,Reg)=[R.SSst(:,Reg);R.SSd(:,Reg)];

for i=102:115
[SS,EP]=RegGae_steadystate_crisis(Gamma(:,i));
Gamma(38:42,i)=EP;
[R.T(:,:,i), R.S(:,:,i), R.Q(:,:,i), R.Z(:,:,i), R.C(:,:,i), R.D(:,:,i), ...
    R.SSd(:,i), R.st(i), R.p(i), R.j(i), R.Ast(:,:,:,i), R.SSst(:,i), R.Dst(:,:,i) ]=...
    RegGae_jacobian('SKANK_crisis',SS,Gamma(:,i));
R.SS(:,i)=[R.SSst(:,i); R.SSd(:,i)];
end
clear EP; clear SS

[cellstr(M_.endo_names),array2table(R.SS(:,[10 115]))]
[cellstr(M_.param_names),array2table(Gamma(:,101: 115))]
table(cellstr(M_.endo_names),R.SS(:,10),R.SS(:,101));

%% IV Transition matrices TM: 
for i=1:T;
    iota=Iota(T-i+1);
    Reg=iota; 
if iota==10|iota==115
    [TM.J(:,:,iota), TM.P(:,:,iota), TM.U(:,:,iota), TM.M(:,:,iota), TM.H(:,:,iota), TM.K(:,:,iota), ...
    TM.V(:,:,iota), TM.Y(:,:,iota), TM.W(:,:,iota),TM.SS_iota(:,:,iota) ]=...
        RegGae_infinite(R.T(:,:,Reg), R.S(:,:,Reg), R.Q(:,:,Reg), R.Z(:,:,Reg), R.C(:,:,Reg), R.D(:,:,Reg), ...
            R.SSd(:,Reg), R.st(Reg), R.p(Reg), R.j(Reg), R.Ast(:,:,:,Reg), R.SSst(:,Reg), R.Dst(:,:,Reg));
else 
    [TM.J(:,:,iota), TM.P(:,:,iota), TM.U(:,:,iota), TM.M(:,:,iota), TM.H(:,:,iota), TM.K(:,:,iota), ...
    TM.V(:,:,iota), TM.Y(:,:,iota), TM.W(:,:,iota), TM.SS_iota(:,:,iota),BK]=...
        RegGae_finite(R.T(:,:,Reg), R.S(:,:,Reg), R.Q(:,:,Reg), R.Z(:,:,Reg), R.C(:,:,Reg), R.D(:,:,Reg), ...
            R.SSd(:,Reg), R.st(Reg), R.p(Reg), R.j(Reg), R.Ast(:,:,:,Reg), R.SSst(:,Reg), R.Dst(:,:,Reg), TM.V(:,:,Iota(T-i+1+1)), TM.Y(:,:,Iota(T-i+1+1)));
end;
end;
%

%}
% Miscellaneous
%{
% some controls (just to inspect the non-homogenous, level term in each iota)
[TM.P(:,:,1) TM.P(:,:,2) TM.P(:,:,3) TM.P(:,:,4) TM.P(:,:,5) ;...
    TM.H(:,:,1) TM.H(:,:,2) TM.H(:,:,3) TM.H(:,:,4) TM.H(:,:,5) ; ...
    TM.Y(:,:,1) TM.Y(:,:,2) TM.Y(:,:,3) TM.Y(:,:,4) TM.Y(:,:,5) ] 

SS=[R.SSst(:,:,5) R.SSst(:,:,4) R.SSst(:,:,3) R.SSst(:,:,2) R.SSst(:,:,1); ...
    R.SSd(:,:,5) R.SSd(:,:,4) R.SSd(:,:,3) R.SSd(:,:,2) R.SSd(:,:,1)] % just to inspect the regime steady states

RsSS=[TM.SS_iota(:,:,5) TM.SS_iota(:,:,4) TM.SS_iota(:,:,3) TM.SS_iota(:,:,2) TM.SS_iota(:,:,1)]; ...
    
table([1:R.p(1)+R.j(1)+R.st(1)]', R.names, Steady_states);
table_ss=[[1:R.p(1)+R.j(1)+R.st(1)]', Steady_states];
%}

% III Simulation
%close all

%% V-Simulation: initial state

x=[];
e=randn(M_.exo_nbr,T); 
e(:,1)=vareps;

% First period: normal history
iota=Iota(1);
x(:,1)=[TM.J(:,:,iota); TM.M(:,:,iota); TM.V(:,:,iota)]*R.SSd(1:R.p(iota),iota)+...
    [TM.P(:,:,iota); TM.H(:,:,iota); TM.Y(:,:,iota)]+...
    +[TM.U(:,:,iota); TM.K(:,:,iota); TM.W(:,:,iota)]*e(:,1);

[R.SS(:,10) x(:,1) [R.SSst(:,10); R.SSd(:,10)]];

iota=10;
%[TM.J(:,:,iota); TM.M(:,:,iota); TM.V(:,:,iota)]-oo_.dr.ghx;

for t=2:T;
    x(:,t)=[TM.J(:,:,Iota(t)); TM.M(:,:,Iota(t)); TM.V(:,:,Iota(t))]*x(R.st(1)+1:R.st(1)+R.p(1),t-1)+...
           [TM.P(:,:,Iota(t)); TM.H(:,:,Iota(t)); TM.Y(:,:,Iota(t))]+...
          [TM.U(:,:,Iota(t)); TM.K(:,:,Iota(t)); TM.W(:,:,Iota(t))]*e(:,t)*0;       
end

M_.endo_names(35,:)='LR    ';
x(35,:)=x(17,:)./(x(13,:)+x(15,:)+x(16,:));

M_.endo_names(36,:)='DF    ';
x(36,:)=x(15,:)+x(16,:);

M_.endo_names(37,:)='DFDFss';
x(37,:)=x(36,:)/(R.SS(15,10)+R.SS(16,10));


%% Generate Figures 
close all
v=1;
for f=1:5
    loc=1; figure;
    while loc<10 & v<=size(x,1)
subplot(3,3,loc)
ytickformat('%.2f')
if v==35|v==36|v==37
    plot(0:T-1,x(v,1:T),'b')
else
plot(0:T-1,x(v,1:T),'b',0:T-1,R.SS(v,Iota),'r:') %,0:T-1,TM.SS_iota(v,Iota),'r'
line([0 T-1], [R.SS(v,Iota(1)) R.SS(v,Iota(1))], 'Color','black','LineStyle',':')
end
title(M_.endo_names(v,:))
%legend('location','southoutside')
%legend(M_.endo_names(v,:),'SS', 'SS_{\iota}')

grid on
loc=loc+1; v=v+1;
    end   
end

% saveas(f4,'Figure_04.png')
    
figure;
v=5;
subplot(2,2,1)
plot(0:T-1,x(v,1:T),'b',0:T-1,R.SS(v,Iota),'r:') %,0:T-1,TM.SS_iota(v,Iota),'r'
title(M_.endo_names(v,:))
%legend('location','southoutside')
%legend(M_.endo_names(v,:),'SS', 'SS_{\iota}')
line([0 T-1], [R.SS(v,Iota(1)) R.SS(v,Iota(1))], 'Color','black','LineStyle',':')
grid on

v=26;
subplot(2,2,2)
plot(0:T-1,x(v,1:T),'b',0:T-1,R.SS(v,Iota),'r:') % ,0:T-1,TM.SS_iota(v,Iota),'r'
title(M_.endo_names(v,:))
%legend('location','southoutside')
%legend(M_.endo_names(v,:),'SS', 'SS_{\iota}')
line([0 T-1], [R.SS(v,Iota(1)) R.SS(v,Iota(1))], 'Color','black','LineStyle',':')
grid on

%}
asset=[x(16,:);x(15,:);x(13,:)];
liab=[-x(17,:);-x(14,:);-Gamma(7,Iota).*x(23,:)];

subplot(2,2,3)
title('Balance sheet')
bar(1:T,asset','stacked')
hold on
bar(1:T,liab','stacked')
hold off