function [eta1,eta2,nu1,nu2] = updateGlobals(alpha, eta1, rho1, eta2, rho2, nu1, rho_nu1, nu2, rho_nu2, x, beta, N_t, lambda_t_s)
    eta1 = eta1 + alpha*(1-sum(x,2)+rho1);
    eta1(eta1<0)=0;
    eta2 = eta2 + alpha*(sum(x,2)-N_t+rho2);
    eta2(eta2<0)=0;
    nu1 = nu1 + alpha*(lambda_t_s-sum(beta,3)+rho_nu1 - 0.3);
    nu1(nu1<0)=0;
    nu2 = nu2 + alpha*(sum(beta,3)-lambda_t_s+rho_nu2 + 0.3);
    nu2(nu2<0)=0;
end

