clc
clear all
close all
tic
numNode = 200;
numEdge = numNode-1;
numVariable = 3;
Kg = zeros(1,numNode-1);
for  branchIndex = 1: numNode - 1
    el=canonicalNets(numNode,'tree',branchIndex);
    adjMatrix = edgeL2adj(el);
    %A1 = zeros(2*numEdge*numVariable,numNode*numVariable);
    A1Hat = zeros(2*numEdge,numNode);
    %A2 = zeros(2*numEdge*numVariable,numNode*numVariable);
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
    %     A1 = kron(A1Hat,eye(numVariable));
    %     A2 = kron(A2Hat,eye(numVariable));
    A1 = A1Hat;
    A2 = A2Hat;
    A = [A1;A2];
    M_p = A1' + A2';
    M_m = A1' - A2';
    sigmaMaxMp = max(svd(M_p));
    M_mSVD = svd(M_m);
    sigmaMin = M_mSVD(end-1);
    Kg(branchIndex) = sigmaMaxMp/sigmaMin;
end
toc