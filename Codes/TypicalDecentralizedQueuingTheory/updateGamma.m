function [ Gamma ] = updateGamma(lambda_t_u_l, K1_t, K2_t, pri_l )
Gamma = zeros(size(lambda_t_u_l));
for i=1:size(lambda_t_u_l,3) %l_l+1
    Gamma(:,:,i) = pri_l(i)*(K1_t.*lambda_t_u_l(:,:,i)+K2_t);
end

end

