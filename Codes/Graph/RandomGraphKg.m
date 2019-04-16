clc
clear all
close all
tic
numNode = 100;
connectivityRatio = .1;
numEdge = round(connectivityRatio*numNode*(numNode-1)/2);
numVariable = 3;
numConnected = 0;
numExpriment = 100;
Kg = zeros(1,numExpriment);
for  exp = 1: numExpriment
    %el=randomGraph(numNode,.2,numEdge);
    %adjMatrix = edgeL2adj(el);
    adjMatrix = randomGraph(numNode,.2,numEdge);
    flag = isConnected(adjMatrix);
    if flag == 1
        numConnected = numConnected + 1;
        A1Hat = zeros(2*numEdge,numNode);
        A2Hat = zeros(2*numEdge,numNode);
        z = zeros(2*numEdge,2);
        kk = 1;
        for ii = 1:numNode
            for jj=1:numNode
                if adjMatrix(ii,jj)==1
                    z(kk,:) = [ii,jj];
                    kk = kk+1;
                end
            end
        end
        for qq = 1:2*numEdge
            A1Hat(qq,z(qq,1)) = 1;
            A2Hat(qq,z(qq,2)) = 1;
        end
        A1 = A1Hat;
        A2 = A2Hat;
        A = [A1;A2];
        M_p = A1' + A2';
        M_m = A1' - A2';
        sigmaMaxMp = max(svd(M_p));
        M_mSVD = svd(M_m);
        sigmaMin = M_mSVD(end-1);
        Kg(numConnected) = sigmaMaxMp;
        if Kg(numConnected)== 0
            break
        end
    end
end
connectedPercent = numConnected/numExpriment
if numConnected > 0
    hist(Kg(1:numConnected))
end
% fileID = ['Numedge',num2str(numEdge),'NumNode',num2str(numNode)]
% save(fileID)
toc