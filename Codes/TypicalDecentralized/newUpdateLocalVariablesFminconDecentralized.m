function [ x ] = newUpdateLocalVariablesFminconDecentralized( Cost_m, lambda, nu, tau_m, delta, l_l, l_t, w_t, cap_m)

f = @(x) x*Cost_m + x*diag(lambda*tau_m') + x*nu - (lambda'*delta + sum(nu))/3/l_l;
x = fmincon(f, zeros(1,l_t), w_t', cap_m, [], [], zeros(1,l_t), ones(1,l_t)); 

end