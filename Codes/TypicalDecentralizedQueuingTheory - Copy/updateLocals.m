function [x, beta, gamma] = updateLocals(eta1, eta2, nu, topo, cap, tauTr, pri)
    l_t = topo.l_t;
    l_s = topo.l_s;
    lambda_t_s = topo.lambda_t_s;
    epsilon = topo.epsilon;
    Q_t = topo.Q_t;
    k1_t = topo.k1_t;
    k2_t = topo.k2_t;
    w_t = topo.w_t;
    CPU = topo.CPU;
    delta_t = topo.delta_t;
    K2_t = sum(k2_t,2);
    K1_t = sum(k1_t,2);
    
    PSI_RENG = 1 : l_t;
    PSI_RENG_MAT = reshape(PSI_RENG, l_t, 1);
    X_RENG = l_t+1 : 2*l_t;
    X_RENG_MAT = reshape(X_RENG, l_t, 1);
    GAMMA_RENG = 2*l_t+1 : (2+l_s)*l_t;
    GAMMA_RENG_MAT = reshape(GAMMA_RENG, l_t, l_s);
    BETA_RENG = (2+l_s)*l_t+1 : (2+2*l_s)*l_t;
    BETA_RENG_MAT = reshape(BETA_RENG, l_t, l_s);
    PHI_RENG = (2+2*l_s)*l_t+1 : (2+3*l_s)*l_t;
    PHI_RENG_MAT = reshape(PHI_RENG, l_t, l_s);
    numofTotVar = (2+3*l_s)*l_t;
    %%
    blx = zeros(1, numofTotVar);
    bux = inf*ones(1, numofTotVar);
    aRow = zeros(1,numofTotVar);
    A = [];
    buc = [];
    blc = [];
    numofCons=0;
    %% constraint 5
    bux(BETA_RENG_MAT) = lambda_t_s;
    bux(X_RENG) = 1;
    bux(GAMMA_RENG) = 1;
    %% constraint 7
    for k=1:size(GAMMA_RENG_MAT,2) %l_s
        for i=1:size(GAMMA_RENG_MAT,1) %l_t
            numofCons = numofCons + 1;
            aRow = zeros(1,numofTotVar);
            aRow(GAMMA_RENG_MAT(i,k)) = 1;
            aRow(X_RENG_MAT(i)) = -1;
            A(numofCons,:) = aRow;
            buc(numofCons) = 0;
            blc(numofCons) = -inf;
        end
    end    
    %% constraint 8
    for k=1:size(BETA_RENG_MAT,2) %l_s
        for i=1:size(BETA_RENG_MAT,1) %l_t
            numofCons = numofCons + 1;
            aRow = zeros(1,numofTotVar);
            aRow(BETA_RENG_MAT(i,k)) = 1;
            aRow(GAMMA_RENG_MAT(i,k)) = -lambda_t_s(i,k);
            A(numofCons,:) = aRow;
            buc(numofCons) = 0;
            blc(numofCons) = -inf;
            
            numofCons = numofCons + 1;
            aRow = zeros(1,numofTotVar);
            aRow(BETA_RENG_MAT(i,k)) = 1;
            aRow(GAMMA_RENG_MAT(i,k)) = -1;
            A(numofCons,:) = aRow;
            buc(numofCons) = inf;
            blc(numofCons) = -1+epsilon;
        end
    end    
    %% constraint 11
    for i=1:size(PSI_RENG_MAT,1) %l_t
        numofCons = numofCons + 1;
        aRow = zeros(1,numofTotVar);
        aRow(PSI_RENG_MAT(i)) = 1;
        aRow(BETA_RENG_MAT(i,:)) = -1;
        A(numofCons,:) = aRow;
        buc(numofCons) = 0;
        blc(numofCons) = -inf;
        
        numofCons = numofCons + 1;
        aRow = zeros(1,numofTotVar);
        aRow(PSI_RENG_MAT(i)) = 1;
        aRow(X_RENG_MAT(i)) = -Q_t(i);
        A(numofCons,:) = aRow;
        buc(numofCons) = 0;
        blc(numofCons) = -inf;
        
        numofCons = numofCons + 1;
        aRow = zeros(1,numofTotVar);
        aRow(PSI_RENG_MAT(i)) = 1;
        aRow(BETA_RENG_MAT(i,:)) = -1;
        aRow(X_RENG_MAT(i)) = -Q_t(i);
        A(numofCons,:) = aRow;
        buc(numofCons) = inf;
        blc(numofCons) = -Q_t(i);
    end

    %% constraint 12
    for r=1:size(k1_t,2) %l_r
        numofCons = numofCons + 1;
        aRow = zeros(1,numofTotVar);
        aRow(PSI_RENG_MAT(:)) = k1_t(:,r);
        aRow(X_RENG_MAT(:)) = k2_t(:,r);
        A(numofCons,:) = aRow;
        buc(numofCons) = cap(r);
        blc(numofCons) = -inf;
    end

    %% constraint 14
    for k=1:size(PHI_RENG_MAT,2) %l_s
        for i=1:size(PHI_RENG_MAT,1) %l_t
            numofCons = numofCons + 1;
            aRow = zeros(1,numofTotVar);
            aRow(PHI_RENG_MAT(i,k)) = (k1_t(i,CPU)-w_t(i))*tauTr;
            aRow(GAMMA_RENG_MAT(i,k)) = k2_t(i,CPU)*tauTr + w_t(i);
            aRow(BETA_RENG_MAT(i,:)) = -(k1_t(i,CPU)-w_t(i))*delta_t(i);
            A(numofCons,:) = aRow;
            buc(numofCons) = k2_t(i,CPU)*delta_t(i);
            blc(numofCons) = -inf;
            
            numofCons = numofCons + 1;
            aRow = zeros(1,numofTotVar);
            aRow(PHI_RENG_MAT(i,k)) = 1;
            aRow(BETA_RENG_MAT(i,:)) = -1;
            A(numofCons,:) = aRow;
            buc(numofCons) = 0;
            blc(numofCons) = -inf;
            
            numofCons = numofCons + 1;
            aRow = zeros(1,numofTotVar);
            aRow(PHI_RENG_MAT(i,k)) = 1;
            aRow(GAMMA_RENG_MAT(i,k)) = -Q_t(i);
            A(numofCons,:) = aRow;
            buc(numofCons) = 0;
            blc(numofCons) = -inf;
        
            numofCons = numofCons + 1;
            aRow = zeros(1,numofTotVar);
            aRow(PHI_RENG_MAT(i,k)) = 1;
            aRow(BETA_RENG_MAT(i,:)) = -1;
            aRow(GAMMA_RENG_MAT(i,k)) = -Q_t(i);
            A(numofCons,:) = aRow;
            buc(numofCons) = inf;
            blc(numofCons) = -Q_t(i);
            
        end
    end

    %% constraint 16
    for i=1:size(X_RENG_MAT,1) %l_t
        numofCons = numofCons + 1;
        aRow = zeros(1,numofTotVar);
        aRow(X_RENG_MAT(i)) = epsilon;
        aRow(BETA_RENG_MAT(i,:)) = w_t(i)-k1_t(i);
        A(numofCons,:) = aRow;
        buc(numofCons) = k2_t(i);
        blc(numofCons) = -inf;
    end

    %% Objective
    c = zeros(1,numofTotVar);
    c(BETA_RENG_MAT) = nu;
    c(X_RENG_MAT) = pri*K2_t-eta1+eta2;
    c(PSI_RENG_MAT) = pri*K1_t;
    %% integer variables
    ints = [X_RENG, GAMMA_RENG];
    %% Problem Defenition
    clear prob
    prob.c = c;
    prob.a = A;
    prob.blc = blc;
    prob.buc = buc;
    prob.blx = blx;
    prob.bux = bux;
    prob.ints.sub = ints;
    %% Problem Solving
    % Optimize the problem.
    [~,res] = mosekopt('minimize',prob);

    try 
      % Get the optimal solution.
      x = res.sol.int.xx(X_RENG_MAT);
%       x = res.sol.bas.xx(X_RENG_MAT);
      beta = res.sol.int.xx(BETA_RENG_MAT);
%       beta = res.sol.bas.xx(BETA_RENG_MAT);
      gamma = res.sol.int.xx(GAMMA_RENG_MAT);
%       gamma = res.sol.bas.xx(GAMMA_RENG_MAT);
    catch
      fprintf('MSKERROR: Could not get solution')
    end

end