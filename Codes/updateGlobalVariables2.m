function y = updateGlobalVariables2(alpha, nu, x_e, x_f, x_c)
y = nu + alpha*(sum(x_e)' + sum(x_f)' + sum(x_c)' - 1);
end

