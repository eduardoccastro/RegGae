function [T, S, Q, Z, C, Dd, SSd, st, p, j, Ast, SSst, Dst]=RegGae_jacobian(model,SS,param_values)
%% RegGae: a toolkit for macroprudencial policy
% Copyright 2020, Eduardo C. Castro ecastro@bcb.gov.br
% May be used at your own risk but proper credit is required
% This function: runs Dynare for 'model', extracts its Jacobian, steady-state and shock matrices.

% Output notation:
% QTZ: "forward" matrix
% QSZ: "current" matrix
% C: non-homogenous term 
% D: matrix premultiplying shocks
% SSd: steady-states of dynamic vars 
% p: number of predetermined variables
% j: number of jumpers (forward-looking variables)
% st: number of static vars
% Ast: matrices for static variables
% SSst: steady-states of static vars
% Dst: shocks matrix for static variables


global oo_ M_

L=M_.lead_lag_incidence';

% Step 0: Construct vector of steady-state values in the order of M_.lead_lag_incidence
% (to recoup the Jacobian below)

steady_state_values=[];
for i=1:3
    for j=1:size(L,1)
        if L(j,i)>0
steady_state_values(L(j,i),1)=SS(j); 
        else 
        end 
    end
end 

% Check, with L, that variable ordering is: static, followed by
% predetermined and lastly jumpers. This was ensured at the dynared mod.
% file.

% Collects the Jacobian matrix g1 from dynare in modfile_dynamic.m
% evaluated at the steady state.
x=zeros(M_.exo_nbr); 

%
[~, g1] =feval(strcat(model,'_dynamic'), steady_state_values, x, param_values, SS, 1); 
%
p=M_.npred; % includes lead-lag variables?
j=M_.nsfwrd; % only pure jumpers
st=M_.nstatic;
n=p+j+st;

% Converting g1 into ABCD format: A.E_t[x^p_t x^j_{t+1}]'=B[x^p_{t-1} x^j_t]'+C+D.e_t
% where x^p_t are the predetermined vars and x^j_t are the jumpers (mixed variables
% were replaced with auxiliary vars).

% Step 1: Obtaining the Jacobian matrices of the model from g1
% Assembling array A with 3 pages: A(1)E_tx_{t+1}+A(2)x_t+A(3)x_{t-1}+SS+De.e_t

A=[];
for i=1:3
for k=1:size(L,1)
    if L(k,i)==0
        A(:,k,4-i)=zeros(size(L,1),1);
    else A(:,k,4-i)=g1(:,L(k,i));
    end 
end 
end

% Creating the matrix De that premultiplies shock vector
De=g1(:,max(max(L))+1:size(g1,2));

% % Step 2: Creating auxiliary variables for mixed variables (with leads and lags):
% 
% % Substep 2.1: Identifying which aux vars are needed
 
 aux_vars=[];
 n_aux_vars=0;

% for i=1:size(A,2);
%     if sum(A(:,i,1) ~= zeros(n,1))>0 & sum(A(:,i,3) ~= zeros(n,1))>0; % Test: if var exist in t+1 and t-1
%         aux_vars(length(aux_vars)+1)=i;
%         n_aux_vars=n_aux_vars+1;
%     else
%     end
% end
% 
% % Substep: 2.2: Inserting new rows and columns with zeros
% 
% A_3_temp=A(:,:,3); %temp copy of lagged matrix
% SS_temp=SS; %temp copy of SS vector
% A=[A(:,1:p,:), zeros(n, n_aux_vars,3), A(:,p+1:p+j+st,:)]; %middle columns
% SS=[SS(1:p);zeros(n_aux_vars,1);SS(p+1:p+j+st)]; % middle rows of SS vector
% A=[A(:,:,:);zeros(n_aux_vars, size(A,2),3)]; % bottom rows
% 
% % Substep 2.3: Entering new equations' coefficients
% for i=1:length(aux_vars);
%     A(p+j+st+i,p+i,1)=1; %coeffs in future matrix
%     SS(p+i)=SS_temp(aux_vars(i)); % SSs in SS' vector
%     if aux_vars(i) <= p; %coeffs in current matrix 
%         A(p+j+st+i, aux_vars(i), 2)=1;
%         else A(p+j+st+i, aux_vars(i)+n_aux_vars,2)=1;
%     end
% 
% end
% clear SS_temp;
% 
% % Step 2.4: Adjusting previous equations
% for i=1:length(aux_vars);
%     A(1:p+j+st,p+i,2)=A_3_temp(:,aux_vars(i)); 
% end
% 
% clear A_3_temp;
% p=p+n_aux_vars; % number of predets is added by number of aux vars
% n=p+j+st;

% At this stage, we have:
% 1) array A of 3 pages (lead, current and lag) matrices with STATIC vars too;
% 2) matrix De that premultiplies shock vector;
% 3) SS vector (full) of steady-state values

% Eliminating static variables with QR-decomposition

St=A(:,1:st,2);
[Qs, ~]=qr(St); %Does QR decompositon
QsA(:,:,1)=Qs'*A(:,:,1);
QsA(:,:,2)=Qs'*A(:,:,2);
QsA(:,:,3)=Qs'*A(:,:,3);
Ast=QsA(1:st,:,:);

A_lead=Qs'*A(:,:,1);
A_plus=A_lead(:,p+st+1:n);
A_plus_tilde=A_plus(st+1:n,:);

A_lag=Qs'*A(:,:,3);
A_minus=A_lag(:,st+1:st+p);
A_minus_tilde=A_minus(st+1:n,:);

A_0=Qs'*A(:,:,2);
A_0_tilde_plus=A_0(st+1:n,st+p+1:n);
A_0_tilde_minus=A_0(st+1:n,st+1:st+p);

QsDe=Qs'*De;
Dst=QsDe(1:st,:); % only the static variables
Dd=QsDe(st+1:n,:); %only the dynamic variables

% QZ factorization
D=[-A_0_tilde_minus -A_plus_tilde]; 
E=[A_minus_tilde A_0_tilde_plus]; 

C=(D-E)*SS(st+1:n);

[S,T,QQ,ZZ] = qz(E, D,'real');
[S,T,QQ,ZZ] = ordqz(S,T,QQ,ZZ, 'udi') ;
BK_QZ=ordeig(S,T);

Z=ZZ';
Q=QQ';

SSd=SS(st+1:n,:); % only the dynamic steady states.
SSst=SS(1:st,:); % only the static steady states.

end