function [ x ] = newUpdateLocalVariablesFminconDecentralized( C, lambda, nu, tau, delta, l_n, l_t, w, c)

f = @(x) x*C + x*diag(lambda*tau') + x*nu - (lambda'*delta + sum(nu))/3/l_n;
x = fmincon(f, zeros(1,l_t), w', c, [], [], zeros(1,l_t), ones(1,l_t)); 

end

