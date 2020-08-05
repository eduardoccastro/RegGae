function [J, P, U, M, H, K, V, Y, W, SS_iota] = RegGae_infinite(T, S, Q, Z, C, D, SSd, st, p, j, Ast, SSst, Dst)
%% RegGae: a toolkit for macroprudencial policy
% Copyright 2020, Eduardo C. Castro ecastro@bcb.gov.br
% May be used at your own risk but proper credit is required
% This function: computes transition matrices for infinite regimes

nd=p+j; %number of dynamic variables
n=nd+st;

Z22=Z(p+1:nd,p+1:nd);
Z21=Z(p+1:nd,1:p);
S22=S(p+1:nd,p+1:nd);
T22=T(p+1:nd,p+1:nd);

C_til=Q'*C;
C_til_p=C_til(1:p);
C_til_j=C_til(p+1:nd);

D_til=Q'*D;
D_til_p=D_til(1:p,:);
D_til_j=D_til(p+1:nd,:);

V=-inv(Z22)*Z21;
%Y=-inv(Z22)*inv(eye(j)-inv(S22)*T22)*inv(S22)*C_til_j;
Y=-V*SSd(1:p)+SSd(p+1:nd);
W=-inv(Z22)*inv(S22)*D_til_j;

BK_infinite=eig(inv(S22)*T22);

TpZ=T(1:p,:)*Z;
TpZ1=TpZ(:,1:p);
TpZ2=TpZ(:,1+p:nd);
SpZ=S(1:p,:)*Z;
SpZ1=SpZ(:,1:p);
SpZ2=SpZ(:,p+1:nd);

M=inv(TpZ1+TpZ2*V)*(SpZ1+SpZ2*V);
%H=inv(TpZ1+TpZ2*V)*((SpZ2-TpZ2)*Y+C_til_p);
H=-M*SSd(1:p,:)+SSd(1:p);
K=inv(TpZ1+TpZ2*V)*(SpZ2*W+D_til_p);

%% Recovering static variables

A_plus=Ast(:,st+p+1:n,1);
A_0s=Ast(:,1:st,2);
A_0_minus=Ast(:,st+1:st+p,2);
A_0_plus=Ast(:,st+p+1:n,2);
A_minus=Ast(:,st+1:st+p,3);

J=-inv(A_0s)*(A_plus*V*M+A_0_minus*M+A_0_plus*V+A_minus);
P=SSst-J*SSd(1:p,:);
U=-inv(A_0s)*(A_plus*V*K+A_0_minus*K+A_0_plus*W+Dst);

%% Proxy steady states (iota)

xp_iota=inv(eye(p)-M)*H;
xj_iota=V*xp_iota+Y;
xst_iota=J*xp_iota+P;
SS_iota=[xst_iota; xp_iota;xj_iota];

end

