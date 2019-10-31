function [ c, ceq ] = consDecentralized(x, m, cap_m, w_t, tau_m ,delta_t)

c = [x(:,m)'*w_t-cap_m, (sum((x.*tau_m'),2)-delta_t)'];
ceq = sum(x,2) - 1;

end

