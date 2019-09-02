function [ x ] = updateLocalVariablesFminconDecentralized( m, Cost_m, v_m, z, rho, l_t, l_m, cap_m, w_t, tau_m, delta_t) 

f = @(x) x(:,m)'*Cost_m + sum(sum(v_m.*(x-z))) + sum(sum((x-z).*(x-z)))*rho/2;

cons = @(x) consDecentralized(x, m, cap_m, w_t, tau_m ,delta_t);

x = fmincon(f, zeros(l_t,l_m), [],[],[],[],zeros(l_t,l_m), ones(l_t,l_m), cons);
end

