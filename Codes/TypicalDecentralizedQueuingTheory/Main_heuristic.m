addpath D:\Installed\Mosek\9.1\toolbox\r2015a
clear
clc

l_c = 2;
l_f = 3;
l_e = 4;
l_s = 3;
l_t = 6;%%%%%increase l_t
l_l = l_c + l_f + l_e;

N_r = 3; %index1:CPU, index2:Storage, index3:RAM  
CPU = 1;
STRG = 2;
RAM = 3;

cap_c = 100+rand(l_c, N_r);
cap_f = 90+rand(l_f, N_r);
cap_e = 80+rand(l_e, N_r);
cap_l = [cap_c; cap_f; cap_e];

delta_t = 30+rand(l_t,1);%%%%
w_t = 4+rand(l_t,1);%%%%%%%%
N_t = 3+zeros(l_t,1); %%%%more better

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

pri_c = 0.4+rand(l_c, 1);
pri_f = 0.4+rand(l_f, 1);
pri_e = 0.9+rand(l_e, 1);
pri_l = [pri_c; pri_f; pri_e];

tauTr_c = 10+rand(l_c, 1);
tauTr_f = 4+rand(l_f, 1);
tauTr_e = 2+rand(l_e, 1);
tauTr_l = [tauTr_c; tauTr_f; tauTr_e];

lambda_t_s = 5+rand(l_t, l_s);%%%%%%%%
epsilon = 0.001;
Q_t = sum(lambda_t_s, 2);
topo = struct('l_t', l_t, 'l_s', l_s, 'lambda_t_s', lambda_t_s, 'epsilon', epsilon, 'Q_t', Q_t, 'k1_t', k1_t, 'k2_t', k2_t, 'w_t', w_t, 'CPU', CPU, 'delta_t', delta_t);

L_v = sum(N_t);
%% Heuristic
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

M_t = 10000;

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
            Gamma_t_u_l = updateGamma( lambda_t_u_l, K1_t, K2_t, [0;pri_l] );
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
                        if j == 0
                            Theta_n_i_j(XP_n_j==i) = phi_n_prev(XP_n_j==i);
                            if i == 0
                                T_n_i_j(XP_n_j==i) = 0;
                            else
                                T_n_i_j(XP_n_j==i) = 0;%%%%%%%%%
                            end
                        else           
                            Theta_n_i_j(XP_n_j==i) = Gamma_t_u_l(t,u,j+1) + phi_n_prev(XP_n_j==i);
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

Gamma_t_u_l = updateGamma( lambda_t_u_l, K1_t, K2_t, [0;pri_l] );

x = zeros(l_t, l_l);
for n = 1 : H
    t = tt(n);
    u = uu(n);
    
    if P(n) ~= 0
        x(t,P(n)) = 1;
    end
    
    if tau_t_u_l(t,u,P(n)+1) > delta_t(t)
        if P(n) ~= 0
            display('jjjjjjjjjjjjj')
            x(t,P(n)) = 0;
        end
    end    
end

for t = 1:l_t
    for l = 1 : l_l
        x(t,u)*Gamma_t_u_l(t,
    end
end
P
H
lambda_t_u_l(:,:,1)
x