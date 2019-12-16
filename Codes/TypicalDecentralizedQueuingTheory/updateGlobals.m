function [eta1,eta2,nu] = updateGlobals(alpha, eta1, eta2, nu, x, beta, N_t, lambda_t_s)
    eta1 = eta1 + alpha*(1-sum(x,2));
    eta1(eta1<0)=0;
    eta2 = eta2 + alpha*(sum(x,2)-N_t);
    eta2(eta2<0)=0;
    nu = nu + alpha*(lambda_t_s-sum(beta,3));
end

