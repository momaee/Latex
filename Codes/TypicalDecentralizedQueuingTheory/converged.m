addpath(genpath('Graph'))
addpath D:\Installed\Mosek\9.1\toolbox\r2015a
clear
clc

l_c = 2;
l_f = 3;
l_e = 4;
l_s = 4;
l_t = 4;
l_n = l_c + l_f + l_e;

N_r = 3; %index1:CPU, index2:Storage, index3: RAM  
CPU = 1;
STRG = 2;
RAM = 3;

cap_c = 112+rand(l_c, N_r);
cap_f = 114+rand(l_f, N_r);
cap_e = 116+rand(l_e, N_r);
cap_n = [cap_c; cap_f; cap_e];

delta_t = 12+zeros(l_t,1);%%%%
w_t = 1+rand(l_t,1);%%%%%%%%
N_t = 2+rand(l_t,1);

k1_t = zeros(l_t,N_r);
k2_t = zeros(l_t,N_r);
k1_t(:,CPU) = w_t + 1;
k2_t(:,CPU) = k1_t(:,CPU) .* rand(l_t,1);
k1_t(:,RAM) = 5*rand(l_t,1);
k2_t(:,RAM) = 5*rand(l_t,1);
K1_t = sum(k1_t,2);
K2_t = sum(k2_t,2);

pri_c = 1+rand(l_c, 1);
pri_f = 8+rand(l_f, 1);
pri_e = 13+rand(l_e, 1);
pri_n = [pri_c; pri_f; pri_e];

tauTr_c = 6+rand(l_c, 1);
tauTr_f = 4+rand(l_f, 1);
tauTr_e = 2+rand(l_e, 1);
tauTr_n = [tauTr_c; tauTr_f; tauTr_e];

lambda_t_s = 13+zeros(l_t, l_s);
epsilon = 0.01;