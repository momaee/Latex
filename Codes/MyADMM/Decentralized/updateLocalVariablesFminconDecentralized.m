function [ b, d, u ] = updateLocalVariablesFminconDecentralized( i, y_i, z_k , betaMin, betaMax, deltaMin, deltaMax, Arow_i, rho, DeltaTilda, epsilonBar)

f = @(x) (pow_p(x(1),-1)-betaMax^(-1))/(betaMin^(-1)-betaMax^(-1));
g = @(x) ((x(2) - DeltaTilda)^(-1)-(1-deltaMin)^(-1))/((1-deltaMax)^(-1)-(1-deltaMin)^(-1));
h = @(x) (y_i')*(x(3:end)'-z_k) + rho/2*sum((x(3:end)' - z_k).^2);
l=@(x) f(x)+g(x)+h(x);

c = @(x) consDecentralized(x,i,betaMin, betaMax, deltaMin, deltaMax, Arow_i, DeltaTilda, epsilonBar);

x = fmincon(l,[betaMin, DeltaTilda+1-deltaMax, 0.5*ones(1,length(Arow_i))],[],[],[],[],[],[],c);

b = x(1);
d = DeltaTilda+1-x(2);
u = x(3:end);

end