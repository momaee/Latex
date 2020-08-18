function [ tau ] = updateTau( lambda_t_u_l, k1_t, k2_t, w_t, tauTr_l )

tau = zeros(size(lambda_t_u_l));%%%%%%%%%%index 1
for i=2:size(lambda_t_u_l,3) %l_l+1
    tau(:,:,i) = tauTr_l(i) + w_t./((k1_t-w_t).*lambda_t_u_l(:,:,i)+k2_t);
end

end

