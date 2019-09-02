function [c,ceq] = consDistributed(x,i,betaMin, betaMax,deltaMin, deltaMax, Arow_i, DeltaTilda, epsilonBar)
c = [x(2)-(DeltaTilda+1-deltaMin), (DeltaTilda+1-deltaMax)-x(2), x(1)-betaMax, betaMin-x(1), (x(1)*(Arow_i*x(3:end)')+x(2)*x(i+2))/((DeltaTilda+1-epsilonBar)*x(i+2)) - 1, 0.0001-x(3:end)];
ceq = prod(x(3:end)) - 1;
% ceq = [];