addpath(genpath('Graph'))
addpath D:\Installed\Mosek\9.1\toolbox\r2015a
clear
clc

l_c = 2;
l_f = 3;
l_e = 4;
l_s = 3;
l_t = 6;
l_n = l_c + l_f + l_e;

N_r = 3; %index1:CPU, index2:Storage, index3:RAM  
CPU = 1;
STRG = 2;
RAM = 3;

cap_c = 100+rand(l_c, N_r);
cap_f = 90+rand(l_f, N_r);
cap_e = 80+rand(l_e, N_r);
cap_n = [cap_c; cap_f; cap_e];

delta_t = 30+rand(l_t,1);%%%%
w_t = 4+rand(l_t,1);%%%%%%%%
N_t = 6+zeros(l_t,1); %%%%more better

k1_t = zeros(l_t,N_r);
k2_t = zeros(l_t,N_r);
k1_t(:,CPU) = w_t + 1;
k2_t(:,CPU) = k1_t(:,CPU) .* rand(l_t,1);
k1_t(:,STRG) = 5*rand(l_t,1);
k2_t(:,STRG) = 5*rand(l_t,1);
k1_t(:,RAM) = 5*rand(l_t,1);
k2_t(:,RAM) = 5*rand(l_t,1);
K1_t = sum(k1_t,2);
K2_t = sum(k2_t,2);

pri_c = 0.1+rand(l_c, 1);
pri_f = 0.4+rand(l_f, 1);
pri_e = 0.9+rand(l_e, 1);
pri_n = [pri_c; pri_f; pri_e];

tauTr_c = 10+rand(l_c, 1);
tauTr_f = 4+rand(l_f, 1);
tauTr_e = 2+rand(l_e, 1);
tauTr_n = [tauTr_c; tauTr_f; tauTr_e];

lambda_t_s = 5+rand(l_t, l_s);%%%%%%%%
epsilon = 0.001;
Q_t = sum(lambda_t_s, 2);
topo = struct('l_t', l_t, 'l_s', l_s, 'lambda_t_s', lambda_t_s, 'epsilon', epsilon, 'Q_t', Q_t, 'k1_t', k1_t, 'k2_t', k2_t, 'w_t', w_t, 'CPU', CPU, 'delta_t', delta_t);
%% Centralized
PSI_RENG = 1 : l_t*l_n;
PSI_RENG_MAT = reshape(PSI_RENG, l_t, l_n);
X_RENG = l_t*l_n+1 : 2*l_t*l_n;
X_RENG_MAT = reshape(X_RENG, l_t, l_n);
BETA_RENG = 2*l_t*l_n+1 : (2+l_s)*l_t*l_n;
BETA_RENG_MAT = reshape(BETA_RENG, l_t, l_n, l_s);

numofTotVar = length(PSI_RENG) + length(X_RENG) + length(BETA_RENG);

blx = zeros(1, numofTotVar);
bux = inf*ones(1, numofTotVar);
aRow = zeros(1,numofTotVar);
A = [];
buc = [];
blc = [];
numofCons=0;
%% constraint 5
for k=1:size(lambda_t_s,2) %l_s
    for i=1:size(lambda_t_s,1) %l_t
        bux(BETA_RENG_MAT(i,:,k)) = lambda_t_s(i,k);
    end
end
bux(X_RENG) = 1;

%% constraint 7
for j=1:size(X_RENG_MAT,2) %l_n
    for i=1:size(X_RENG_MAT,1) %l_t
        numofCons = numofCons + 1;
        aRow = zeros(1,numofTotVar);
        aRow(BETA_RENG_MAT(i,j,:)) = 1;
        aRow(X_RENG_MAT(i,j)) = -Q_t(i);
        A(numofCons,:) = aRow;
        buc(numofCons) = 0;
        blc(numofCons) = -inf;
    end
end   

%% constraint 8
for k=1:size(lambda_t_s,2) %l_s
    for i=1:size(lambda_t_s,1) %l_t
        numofCons = numofCons + 1;
        aRow = zeros(1,numofTotVar);
        aRow(BETA_RENG_MAT(i,:,k)) = 1;
        A(numofCons,:) = aRow;
        buc(numofCons) = lambda_t_s(i,k);
        blc(numofCons) = lambda_t_s(i,k);
    end
end
%% constraint 10
for j=1:size(PSI_RENG_MAT,2) %l_n
    for i=1:size(PSI_RENG_MAT,1) %l_t
        numofCons = numofCons + 1;
        aRow = zeros(1,numofTotVar);
        aRow(PSI_RENG_MAT(i,j)) = 1;
        aRow(BETA_RENG_MAT(i,j,:)) = -1;
        A(numofCons,:) = aRow;
        buc(numofCons) = 0;
        blc(numofCons) = -inf;
        
        numofCons = numofCons + 1;
        aRow = zeros(1,numofTotVar);
        aRow(PSI_RENG_MAT(i,j)) = 1;
        aRow(X_RENG_MAT(i,j)) = -Q_t(i);
        A(numofCons,:) = aRow;
        buc(numofCons) = 0;
        blc(numofCons) = -inf;
        
        numofCons = numofCons + 1;
        aRow = zeros(1,numofTotVar);
        aRow(PSI_RENG_MAT(i,j)) = 1;
        aRow(BETA_RENG_MAT(i,j,:)) = -1;
        aRow(X_RENG_MAT(i,j)) = -Q_t(i);
        A(numofCons,:) = aRow;
        buc(numofCons) = inf;
        blc(numofCons) = -Q_t(i);
    end
end

%% constraint 11
for j=1:size(PSI_RENG_MAT,2) %l_n
    for r=1:size(k1_t,2) %l_r
        numofCons = numofCons + 1;
        aRow = zeros(1,numofTotVar);
        aRow(PSI_RENG_MAT(:,j)) = k1_t(:,r);
        aRow(X_RENG_MAT(:,j)) = k2_t(:,r);
        A(numofCons,:) = aRow;
        buc(numofCons) = cap_n(j,r);
        blc(numofCons) = -inf;
    end
end

%% constraint 13
for k=1:size(BETA_RENG_MAT,3) %l_s
    for j=1:size(BETA_RENG_MAT,2) %l_n
        for i=1:size(BETA_RENG_MAT,1) %l_t
            numofCons = numofCons + 1;
            aRow = zeros(1,numofTotVar);
            aRow(PSI_RENG_MAT(i,j)) = (k1_t(i,CPU)-w_t(i))*tauTr_n(j);
            aRow(X_RENG_MAT(i,j)) = k2_t(i,CPU)*tauTr_n(j) + w_t(i);
            aRow(BETA_RENG_MAT(i,j,:)) = -(k1_t(i,CPU)-w_t(i))*delta_t(i);
            A(numofCons,:) = aRow;
            buc(numofCons) = k2_t(i,CPU)*delta_t(i);
            blc(numofCons) = -inf;
        end
    end
end    

%% constraint 14
for i=1:size(X_RENG_MAT,1) %l_t
    numofCons = numofCons + 1;
    aRow = zeros(1,numofTotVar);
    aRow(X_RENG_MAT(i,:)) = 1;
     A(numofCons,:) = aRow;
     buc(numofCons) = N_t(i);
     blc(numofCons) = 1;
end

%% constraint 15
for j=1:size(X_RENG_MAT,2) %l_n
    for i=1:size(X_RENG_MAT,1) %l_t
        numofCons = numofCons + 1;
        aRow = zeros(1,numofTotVar);
        aRow(X_RENG_MAT(i,j)) = epsilon;
        aRow(BETA_RENG_MAT(i,j,:)) = w_t(i)-k1_t(i);
        A(numofCons,:) = aRow;
        buc(numofCons) = k2_t(i);
        blc(numofCons) = -inf;
    end
end
        
        
%% Objective
c = zeros(1,numofTotVar);
c(PSI_RENG_MAT) = K1_t*pri_n';
c(X_RENG_MAT) = K2_t*pri_n';
%% integer variables
ints = [X_RENG];
%% Problem Defenition
clear prob
prob.c = c;
prob.a = A;
prob.blc = blc;
prob.buc = buc;
prob.blx = blx;
prob.bux = bux;
prob.ints.sub = ints;

%% Problem Solving
% Optimize the problem.
[r,res] = mosekopt('minimize',prob);

try 
  % Display the optimal solution.
  res.sol.int.xx(X_RENG_MAT)
catch
  fprintf('MSKERROR: Could not get solution')
end