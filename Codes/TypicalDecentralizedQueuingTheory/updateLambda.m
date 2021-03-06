function [ lambda_t_u_l ] = updateLambda( lambda_t_u_l, t, N, N_t)

if N < N_t
    lambda_t_u_l(t,N+1,:) = lambda_t_u_l(t,N+1,:) + (1/3)*lambda_t_u_l(t,N,:);
    lambda_t_u_l(t,N,:) =  (2/3)*lambda_t_u_l(t,N,:);
end
if N == N_t
    lambda_t_u_l(t,1:N-1,:) = lambda_t_u_l(t,1:N-1,:) + (1/3)*repmat(lambda_t_u_l(t,N,:),1,N-1)/(N-1);
    lambda_t_u_l(t,N,:) =  (2/3)*lambda_t_u_l(t,N,:);
end

end

