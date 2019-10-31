function [adjacencyMat,A1,A2,Mp,Mm,Lp,Lm,W,z] = networkGeneration(treeFalg,branchFactor,regularFlag,degree,numNode,numEdge)
if treeFalg==1
    edgeList = canonicalNets(numNode,'tree',branchFactor);
    adjacencyMat = edgeL2adj(edgeList);
elseif regularFlag==1
    edgeList = kregular(numNode,degree);
    adjacencyMat = edgeL2adj(edgeList);
else
    connectedFlag = 0;
    while connectedFlag==0
        adjacencyMat = randomGraph(numNode,.5,numEdge);
        connectedFlag = isConnected(adjacencyMat);
    end
end
A1 = zeros(2*numEdge,numNode);
A2 = zeros(2*numEdge,numNode);
z = zeros(2*numEdge,2);
kk = 1;
for ii = 1:numNode
    for jj=1:numNode
        if adjacencyMat(ii,jj)==1
            z(kk,:) = [ii,jj];
            kk = kk+1;
        end
    end
end
for qq = 1:2*numEdge
    A1(qq,z(qq,1)) = 1;
    A2(qq,z(qq,2)) = 1;
end
% A1 = kron(A1Hat,eye(numVariable));
% A2 = kron(A2Hat,eye(numVariable));
Mp = A1' + A2';
Mm = A1' - A2';
Lp = .5*(Mp*Mp');
Lm = .5*(Mm*Mm');
W = .5*(Lp + Lm);