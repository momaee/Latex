clear 
clc
l_c = 4;
l_f = 4;
l_e = 4;
l_s = 4;
l_t = 10;
l_m = l_c + l_f + l_e;

Cap_c = (6*rand(l_c, 1));
Cap_f = (4*rand(l_f, 1));
Cap_e = (2*rand(l_e, 1));
Cap_m = [Cap_c;Cap_f;Cap_e];

delta_t = 3*rand(l_t,1);
w_t = 2*rand(l_t,1);

Pri_c = 3*rand(l_c, 1);
Pri_f = 2*rand(l_f, 1);
Pri_e = rand(l_e,1);
Pri_m = [Pri_c;Pri_f;Pri_e];

TauTr_c = 3*rand(l_c, 1);
TauTr_f = 2*rand(l_f, 1);
TauTr_e = rand(l_e,1);
TauTr_m = [TauTr_c;TauTr_f;TauTr_e];

R_c = rand(l_c, 1);
R_f = rand(l_f, 1);
R_e = rand(l_e,1);
R_m = [R_c;R_f;R_e];

tau_c = repmat(w_t', [l_c, 1])./repmat(R_c, [1,l_t]) + TauTr_c;
tau_f = repmat(w_t', [l_f, 1])./repmat(R_f, [1,l_t]) + TauTr_f;
tau_e = repmat(w_t', [l_e, 1])./repmat(R_e, [1,l_t]) + TauTr_e;
tau_m = [tau_c;tau_f;tau_e];

numofIterations = 10;

z = zeros(l_t,l_m,numofIterations);
u = zeros(l_t,l_m,l_m, numofIterations);
v = zeros(size(u));
rho = 0.5;

for t = 1:numofIterations-1
   for m=1:l_m
       [u(:,:,m,t+1)] = updateLocalVariablesFminconDecentralized(m, Pri_m(m)*w_t, v(:,:,m,t), z(:,:,t), rho, l_t, l_m, Cap_m(m), w_t, tau_m, delta_t);
   end
   z(:,:,t+1) = mean(u(:,:,:,t+1),3) + mean(v(:,:,:,t),3)/rho;
   v(:,:,:,t+1) = v(:,:,:,t) + rho*(u(:,:,:,t+1) - z(:,:,t+1));
end
