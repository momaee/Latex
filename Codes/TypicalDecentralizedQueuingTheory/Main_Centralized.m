addpath(genpath('Graph'))
addpath D:\Installed\Mosek\9.1\toolbox\r2015a
clear
clc

l_c = 1+0;
l_f = 2+0;
l_e = 3+0;
l_s = 2+0;
l_t = 5;
l_n = l_c + l_f + l_e;

N_r = 3; %index1:CPU, index2:Storage, index3:RAM  
CPU = 1;
STRG = 2;
RAM = 3;

% cap_c = 100+zeros(l_c, N_r);
% cap_f = 90+zeros(l_f, N_r);
% cap_e = 80+zeros(l_e, N_r);

cap_c = 100+rand(l_c, N_r);
cap_f = 90+rand(l_f, N_r);
cap_e = 80+rand(l_e, N_r);

cap_n = [cap_c; cap_f; cap_e];

% delta_t = 30+zeros(l_t,1);%%%%
% w_t = 4+zeros(l_t,1);%%%%%%%%

delta_t = 30+rand(l_t,1);%%%%
w_t = 4+rand(l_t,1);%%%%%%%%
N_t = 3+zeros(l_t,1); %%%%more better

k1_t = zeros(l_t,N_r);
k2_t = zeros(l_t,N_r);

% k1_t(:,CPU) = w_t + 1;
% k2_t(:,CPU) = k1_t(:,CPU) .* ones(l_t,1)*0.5;
% k1_t(:,STRG) = 5*ones(l_t,1)*0.7;
% k2_t(:,STRG) = 5*ones(l_t,1)*0.6;
% k1_t(:,RAM) = 5*ones(l_t,1)*0.5;
% k2_t(:,RAM) = 5*ones(l_t,1)*0.4;

k1_t(:,CPU) = w_t + 1;
k2_t(:,CPU) = k1_t(:,CPU) .* rand(l_t,1);
k1_t(:,STRG) = 5*rand(l_t,1);
k2_t(:,STRG) = 5*rand(l_t,1);
k1_t(:,RAM) = 5*rand(l_t,1);
k2_t(:,RAM) = 5*rand(l_t,1);

K1_t = sum(k1_t,2);
K2_t = sum(k2_t,2);

% pri_c = 0.1+zeros(l_c, 1);
% pri_f = 0.4+zeros(l_f, 1);
% pri_e = 0.9+zeros(l_e, 1);

pri_c = 0.1+rand(l_c, 1);
pri_f = 0.4+rand(l_f, 1);
pri_e = 0.9+rand(l_e, 1);

pri_n = [pri_c; pri_f; pri_e];

% tauTr_c = 10+zeros(l_c, 1);
% tauTr_f = 4+zeros(l_f, 1);
% tauTr_e = 2+zeros(l_e, 1);

tauTr_c = 10+rand(l_c, 1);
tauTr_f = 4+rand(l_f, 1);
tauTr_e = 2+rand(l_e, 1);

tauTr_n = [tauTr_c; tauTr_f; tauTr_e];

% lambda_t_s = 5+zeros(l_t, l_s);%%%%%%%%

lambda_t_s = 5+rand(l_t, l_s);%%%%%%%%
epsilon = 0.001;
Q_t = sum(lambda_t_s, 2);
topo = struct('l_t', l_t, 'l_s', l_s, 'lambda_t_s', lambda_t_s, 'epsilon', epsilon, 'Q_t', Q_t, 'k1_t', k1_t, 'k2_t', k2_t, 'w_t', w_t, 'CPU', CPU, 'delta_t', delta_t);
%% Centralized
tic
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

%% constraint 15, queue stability
% for j=1:size(X_RENG_MAT,2) %l_n
%     for i=1:size(X_RENG_MAT,1) %l_t
%         numofCons = numofCons + 1;
%         aRow = zeros(1,numofTotVar);
%         aRow(X_RENG_MAT(i,j)) = epsilon*w_t(i);
%         aRow(PSI_RENG_MAT(i,j)) = w_t(i);
%         aRow(BETA_RENG_MAT(i,j,:)) = -k1_t(i,CPU);
%         A(numofCons,:) = aRow;
%         buc(numofCons) = k2_t(i,CPU);
%         blc(numofCons) = -inf;
%     end
% end

for j=1:size(X_RENG_MAT,2) %l_n
    for i=1:size(X_RENG_MAT,1) %l_t
        numofCons = numofCons + 1;
        aRow = zeros(1,numofTotVar);
        aRow(X_RENG_MAT(i,j)) = epsilon*w_t(i)-k2_t(i,CPU);
        aRow(PSI_RENG_MAT(i,j)) = w_t(i)-k1_t(i,CPU);
        A(numofCons,:) = aRow;
        buc(numofCons) = 0;
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
  x_cent = res.sol.int.xx(X_RENG_MAT)
  psi = res.sol.int.xx(PSI_RENG_MAT);
  beta = res.sol.int.xx(BETA_RENG_MAT);
  lambda_cent = sum(beta,3)
  res.sol.int.pobjval
catch
  fprintf('MSKERROR: Could not get solution')
end
toc
% cost = 0;
% cost2 = 0;
% for i = 1 : size(x,1)
%     for j = 1 : size(x,2)
%     cost = cost + x(i,j)*pri_n(j)*K2_t(i) + psi(i,j)*pri_n(j)*K1_t(i);
%     cost2 = cost2 + x(i,j)*pri_n(j)*(K1_t(i)*lambda(i,j)+K2_t(i));
%     end
% end
% cost
% cost2

%% heuristic
l_l = l_n;
cap_l = cap_n;
pri_l = pri_n;
tauTr_l = tauTr_n;

L_v = sum(N_t);
%% Heuristic
tic
Z_n_prev = [];
phi_n_prev = [];
Lambda_n_prev = zeros(l_l+1, L_v);
Lambda_n = Lambda_n_prev;
NodeResourse_n_prev = zeros([size(cap_l), l_l+1]) + cap_l;
ResourceIndicator = true;
RepeatIndicator = false;
H = L_v;

lambda_t_u_l = zeros(l_t, max(N_t), l_l+1);
lambda_total = sum(lambda_t_s,2);
lambda_t_u_l(:,1,:) = reshape(( repmat(lambda_total,l_l+1,1)) ,l_t,1,l_l+1);
lambda_t_u_l_base = lambda_t_u_l;

% Gamma_t_u_l = updateGamma( lambda_t_u_l, K1_t, K2_t, [0;pri_l] );
% tau_t_u_l = updateTau( lambda_t_u_l, k1_t(:,CPU), k2_t(:,CPU), w_t, [0;tauTr_l]);

M_t = 100000;

numof_division = zeros(1, L_v);
for t = 1 : l_t
    while true
        for u = 1 : N_t(t)    
            n = sum(N_t(1:t)) - N_t(t) + u;
            RepeatIndicator = false;
            if u == 1
                %backup for reset
                Z_n_prev_bak = Z_n_prev;
                phi_n_prev_bak = phi_n_prev;
                Lambda_n_prev_bak = Lambda_n_prev;
                NodeResourse_n_prev_bak = NodeResourse_n_prev;
            end
            if u == 1
                Z_n = 1:l_l; 
            else
                Z_n = 0:l_l;
            end
            RemovedStates = zeros(1,length(Z_n)) - 1 ;
            phi_n = zeros(1,length(Z_n));
            NodeResourse_n = NodeResourse_n_prev;
            Gamma_t_u_l = updateGamma( lambda_t_u_l, K1_t, K2_t, [0;pri_l], M_t );
            tau_t_u_l = updateTau( lambda_t_u_l, k1_t(:,CPU), k2_t(:,CPU), w_t, [0;tauTr_l]);
            for j = Z_n 
                OD_n_j = 0;
                XP_n_j = Z_n_prev;
                if j ~= 0
                    for i = XP_n_j
                        if NodeResourse_n_prev(j,:,i+1) > k1_t(t,:)*lambda_t_u_l(t,u,j+1) + k2_t(t,:)
                            %node can process task
                            %I used i+1 in argument because i can be zeros
                        else
                            XP_n_j(XP_n_j==i) = []; %% remove i from XP_n_j OK
                        end
                    end
                    OD_n_j = delta_t(t);
                end
                if ~isempty(XP_n_j)
                    T_n_i_j = zeros(1,length(XP_n_j));
                    Theta_n_i_j = zeros(1,length(XP_n_j));
                    for i = XP_n_j
                        % calculate T_n_i_j
                        % calculate Theta_n_i_j
                        Theta_n_i_j(XP_n_j==i) = Gamma_t_u_l(t,u,j+1) + phi_n_prev(XP_n_j==i);
                        if j == 0
%                             Theta_n_i_j(XP_n_j==i) = phi_n_prev(XP_n_j==i);
                            if i == 0
                                T_n_i_j(XP_n_j==i) = 0;
                            else
                                T_n_i_j(XP_n_j==i) = 0;%%%%%%%%%
                            end
                        else           
%                             Theta_n_i_j(XP_n_j==i) = Gamma_t_u_l(t,u,j+1) + phi_n_prev(XP_n_j==i);
                            if i == 0
                                T_n_i_j(XP_n_j==i) = tau_t_u_l(t,u,j+1);%%%%%%%%%%%
                            else
                                T_n_i_j(XP_n_j==i) = tau_t_u_l(t,u,j+1);
                            end
                        end
                    end
                    % calculate I_n_j
                    [~, I_n_j] = min( M_t*(T_n_i_j - OD_n_j).*((T_n_i_j-OD_n_j)>0) + Theta_n_i_j);
                    % calculate phi_n_j
                    phi_n(Z_n==j) = Theta_n_i_j(I_n_j);
                    % calculate Lambda_n
                    Lambda_n(j+1, n) = j;
                    Lambda_n(j+1, 1:n-1) = Lambda_n_prev(XP_n_j(I_n_j)+1, 1:n-1);
                    % update node resoureces
                    NodeResourse_n(:,:,j+1) = NodeResourse_n_prev(:,:,XP_n_j(I_n_j)+1);
                    if j~= 0 
                        NodeResourse_n(j,:,j+1) = NodeResourse_n(j,:,j+1) - k1_t(t,:)*lambda_t_u_l(t,u,j+1) + k2_t(t,:);
                    end
                else
                    if n ~= 1
                        % add j to RemovedStates OK
                        RemovedStates(Z_n==j) = j;
                    end
                end
                if n == 1
                    if NodeResourse_n(j,:,j+1) > k1_t(t,:)*lambda_t_u_l(t,u,j+1) + k2_t(t,:)
                        phi_n(Z_n==j) = M_t*(tau_t_u_l(t,u,j+1) - OD_n_j).*((tau_t_u_l(t,u,j+1)-OD_n_j)>0) + Gamma_t_u_l(t,u,j+1);
                        Lambda_n(j+1, n) = j; 
                        if j~= 0 
                            NodeResourse_n(j,:,j+1) = NodeResourse_n(j,:,j+1) - k1_t(t,:)*lambda_t_u_l(t,u,j+1) + k2_t(t,:);
                        end
                    else
                       RemovedStates(Z_n==j) = j;
                    end
                end

            end
            if sum(RemovedStates ~= -1)> 0.7*length(Z_n)
                numof_division(n) = numof_division(n)+1;
                if numof_division(n) > 7
                    break;
                end
                RepeatIndicator = true;
                lambda_t_u_l = updateLambda(lambda_t_u_l, t, u, N_t(t));
                %update prevs from backups
                Z_n_prev = Z_n_prev_bak;
                phi_n_prev = phi_n_prev_bak;
                Lambda_n_prev = Lambda_n_prev_bak;
                NodeResourse_n_prev = NodeResourse_n_prev_bak;
                break;
            end
            if RemovedStates == Z_n
                ResourceIndicator = false;
                % calculate H  OK
                H = sum(N_t(1:t)) - N_t(t);
                % Determine Viterbi path using H
                [~, zetta] = min(phi_n);%%%%%%%%
                P = Lambda_n(Z_n(zetta)+1,1:H);    
                break
            else
                %remove all states in RemovedStates from Z_n
                phi_n(Z_n==RemovedStates) = [];
                Z_n(Z_n==RemovedStates) = [];
            end
            %update _prev variables
            Z_n_prev = Z_n;
            phi_n_prev = phi_n;
            Lambda_n_prev = Lambda_n;
            NodeResourse_n_prev = NodeResourse_n;
            
        end %% n
        if RepeatIndicator == false
           break; 
        end
    end
    if ResourceIndicator == false
        break;
    end
end
if ResourceIndicator == true
    % Determine the viterbi path using H
    [~, zetta] = min(phi_n);%%%%%%%%
    P = Lambda_n(Z_n(zetta)+1,1:H);
end

tt = zeros(1,H);
uu = zeros(1,H);
for t = 1 : l_t
    tt( sum(N_t(1:t))-N_t(t)+1:sum(N_t(1:t)) ) = t;
    uu( sum(N_t(1:t))-N_t(t)+1:sum(N_t(1:t)) ) = 1 : N_t(t);
end

tau_t_u_l = updateTau( lambda_t_u_l, k1_t(:,CPU), k2_t(:,CPU), w_t, [0;tauTr_l]);

Gamma_t_u_l = updateGamma( lambda_t_u_l, K1_t, K2_t, [0;pri_l], M_t );

x = zeros(l_t, l_l);
lambda = zeros(l_t,l_l);
for n = 1 : H
    t = tt(n);
    u = uu(n);
    
    if P(n) ~= 0
        x(t,P(n)) = 1;
        lambda(t,P(n)) = lambda(t,P(n)) + lambda_t_u_l(t,u,1);
    end
    
    if tau_t_u_l(t,u,P(n)+1) > delta_t(t)
        if P(n) ~= 0
            disp('jjjjjjjjjjjjj');
            x(t,P(n)) = 0;
            lambda(t,P(n)) = 0;
        end
    end    
end

cost2 = 0;
for i = 1 : size(x,1)
    for j = 1 : size(x,2)
    cost2 = cost2 + x(i,j)*pri_l(j)*(K1_t(i)*lambda(i,j)+K2_t(i));
    end
end

P;
H;
lambda_t_u_l(:,:,1);
x
cost2
toc
