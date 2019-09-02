function [ b, d, u ] = updateLocalVariablesFminconDistributed( i, phi_i, u_k , betaMin, betaMax, deltaMin, deltaMax, Arow_i, rho, DeltaTilda, epsilonBar)
n = length(phi_i);

f = @(x) (pow_p(x(1),-1)-betaMax^(-1))/(betaMin^(-1)-betaMax^(-1));
g = @(x) ((x(2) - DeltaTilda)^(-1)-(1-deltaMin)^(-1))/((1-deltaMax)^(-1)-(1-deltaMin)^(-1));
h = @(x) rho * sum(sum(pow_pos( (repmat(x(3:end)',1,n) - ((u_k(:,i)+u_k)/2) ), 2)).*Arow_i) + (phi_i')*(x(3:end)');
l=@(x) f(x)+g(x)+h(x);

c = @(x) consDistributed(x,i,betaMin, betaMax, deltaMin, deltaMax, Arow_i, DeltaTilda, epsilonBar);

x = fmincon(l,[betaMin, DeltaTilda+1-deltaMin, 0.5*ones(1,length(Arow_i))],[],[],[],[],[],[],c);

b = x(1);
d = DeltaTilda+1-x(2);
u = x(3:end);
end