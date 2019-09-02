function y = updateGlobalVariables(alpha, lambda, x_e, x_f, x_c, delta_t, tau_e, tau_f, tau_c)
y = lambda + alpha*(diag(x_e'*tau_e) + diag(x_f'*tau_f) + diag(x_c'*tau_c) - delta_t);
y(y<0)=0;
end

