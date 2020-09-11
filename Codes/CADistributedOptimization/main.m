%% Graph Generation
addpath(genpath('Distributed'))
addpath(genpath('Decentralized'))
addpath(genpath('..\TypicalDecentralizedQueuingTheory\Graph'))
clear
clc
N = 15;
A = [];
while not (isConnected(A)) 
    A = randomGraph ( N, 0.3 ) ; 
end
G = graph(A);
figure 
plot(G)
title('Network Graph')
nodesForPlot = 1:N;
%% Network Configuration
specRadA = max(eig(A));
deltaMax = 0.8;
deltaMin = 0.08;
epsilonBar = 0.1;
DeltaTilda = max(epsilonBar,deltaMax);
epidTresh = (1-deltaMax)/specRadA;
betaMax = 4*epidTresh;
betaMin = 0.3*betaMax;
max(eig(diag(betaMax*ones(1,N))*A - diag(deltaMin*ones(1,N))))
max(eig(diag(betaMax*ones(1,N))*A - diag(deltaMax*ones(1,N))))
max(eig(diag(betaMin*ones(1,N))*A - diag(deltaMin*ones(1,N))))
max(eig(diag(betaMin*ones(1,N))*A - diag(deltaMax*ones(1,N))))
%% Distributed Algorithm
rho = 0.5;
numofIterations = 30;

phi = zeros(N,N,numofIterations);
phi(:,:,1) = ones(N,N);
betaDist = zeros(N,numofIterations);
deltaDist = zeros(N,numofIterations);
uDist = zeros(N,N,numofIterations);
uDist(:,1) = ones(N,1);
for k = 1:numofIterations-1
    for i = 1:N
        phi(:,i,k+1) = phi(:,i,k) + rho * sum(  A(i,:).*(uDist(:,i,k)-uDist(:,:,k)), 2 ) ; 
        [betaDist(i,k+1), deltaDist(i,k+1), uDist(:,i,k+1)] = updateLocalVariablesFminconDistributed(i, phi(:,i,k+1), uDist(:,:,k), betaMin, betaMax, deltaMin, deltaMax, A(i,:), rho, DeltaTilda, epsilonBar);
    end
end

%% Decentralized Algorithm
y = zeros(N,N,numofIterations);
betaDecent = zeros(N,numofIterations);
deltaDecent = zeros(N,numofIterations);
uDecent = zeros(N,N,numofIterations);
uDecent(:,1) = ones(N,1);
z = zeros(N,numofIterations);
z(:,1) = ones(N,1);
for k = 1:numofIterations-1
    for i = 1:N
        [betaDecent(i,k+1), deltaDecent(i,k+1), uDecent(:,i,k+1)] = updateLocalVariablesFminconDecentralized(i, y(:,i,k), z(:,k), betaMin, betaMax, deltaMin, deltaMax, A(i,:), rho, DeltaTilda, epsilonBar);
    end
    z(:,k+1) = mean(uDecent(:,:,k+1),2) + mean(y(:,:,k),2)/rho;
    y(:,:,k+1) = y(:,:,k) + rho*(uDecent(:,:,k+1)-z(:,k+1)) ; 
end

%% Convergence Factors Calculation
SDist = zeros(numofIterations,1);
SDecent = zeros(numofIterations,1);
phiRecord = zeros(N,numofIterations);
yRecord = zeros(N,numofIterations);
for t = 1:numofIterations
    for i=1:N
        SDist(t) = SDist(t) + sum(sum((uDist(:,i,t)-uDist(:,:,t)).^2) .* A(i,:));
        SDecent(t) = SDecent(t) + sum(sum((uDecent(:,i,t)-uDecent(:,:,t)).^2) .* A(i,:));
        phiRecord(i,t) = norm(phi(:,i,t));
        yRecord(i,t) = norm(y(:,i,t));
    end
end
%% Plotting Results
nodesForPlot = [9,15];
figure
plot(SDist/max(SDist),'-s','LineWidth',2)
% title('Convergence of Distributed and Decentralized Solutions')
xlabel('Iterations','FontSize',12)
ylabel('Normalized error at each iteration','FontSize',12)
hold on
plot(SDecent/max(SDecent),'-*','LineWidth',2)
legend('Distributed', 'Decentralized', 'location', 'best')
grid on
hold off

DFE_ParamDist = max(eig(diag(betaDist(:,end))*A - diag(deltaDist(:,end))));
DFE_ParamDecent = max(eig(diag(betaDecent(:,end))*A - diag(deltaDecent(:,end))));
fprintf('Epsilon is %.4f\n', epsilonBar);
fprintf('DFE parameter of distributed solution that should be less than minus epsilon is: %.4f \n', DFE_ParamDist);
fprintf('DFE parameter of decentralized solution that should be less than minus epsilon is: %.4f \n', DFE_ParamDecent);

figure
grid on
hold on
% title('Convergence of Distributed Solution')
xlabel('Iterations','FontSize',12)
ylabel('Normalized absolute value of local variables','FontSize',12)    
msg = strings(length(nodesForPlot),1);
for i=1:length(nodesForPlot)
    n = nodesForPlot(i);
    plot(10*(phiRecord(n,:)-min(phiRecord(n,:))))
    msg(i) = sprintf('random node %i',n);
end        
legend(msg, 'location', 'best')
hold off

figure 
% title('Convergence of Decentralized Solutions')
xlabel('iterations','FontSize',12)
ylabel('Normalized absolute value of lagrangian coefficients','FontSize',12)
grid on
hold on
msg = strings(length(nodesForPlot),1);
for i=1:length(nodesForPlot)
    n = nodesForPlot(i);
    plot(yRecord(n,:))
    msg(i) = sprintf('random node %i',n);
end
legend(msg, 'location', 'best')
hold off

% figure 
% title('Convergence of Distributed and Decentralized Solutions')
% xlabel('iterations')
% ylabel('beta of nodes')
% hold on
% grid on
% msg = strings(2*length(nodesForPlot),1);
% for i=1:length(nodesForPlot)
%     n = nodesForPlot(i);
%     plot(betaDist(n,:))
%     plot(betaDecent(n,:))
%     msg(2*i-1) = sprintf('random node %i Distributed',n);
%     msg(2*i) = sprintf('random node %i Decentralized',n);
% end
% legend(msg, 'location', 'best')
% hold off

figure 
% title('Convergence of Distributed and Decentralized Solutions')
xlabel('Iterations','FontSize',12)
ylabel('Normalized absolute value of local variables','FontSize',12)    
hold on
grid on
msg = strings(2*length(nodesForPlot),1);
for i=1:length(nodesForPlot)
    n = nodesForPlot(i);
    plot(deltaDist(n,:),'-s','LineWidth',2)
    plot(deltaDecent(n,:),'-*','LineWidth',2)
    msg(2*i-1) = sprintf('random node %i Distributed',n);
    msg(2*i) = sprintf('random node %i Decentralized',n);
end
legend(msg, 'location', 'best')
hold off
