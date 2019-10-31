function [errorNormalized] = ErasureChannelADMM(adjMatrixOrig,rho,numNode,numberVariable,numIteration,A,b,xHatCen)

xHatMatrixOld = zeros(numberVariable,numNode);
xHatMatrixNew = zeros(numberVariable,numNode);
xHatMatrixNew2 = zeros(numberVariable,numNode);
alphaMat = zeros(numberVariable,numNode);
errorNormalized = zeros(numIteration,1);
lastRCVDMat = zeros(numberVariable,numNode,numNode);

for jj = 1:numIteration
    if jj == 1
        for ii = 1:numNode
            nodeIndex = ii;
            neighbor = repmat(adjMatrixOrig(nodeIndex,:),size(xHatMatrixNew,1),1);
            numNeighbor = sum(adjMatrixOrig(nodeIndex,:));
            ACV = A(:,:,nodeIndex);
            bCV = b(:,nodeIndex);
            xHatMatrixNew2(:,nodeIndex) = inv(ACV'*ACV + 2*rho*numNeighbor*eye(numberVariable))*( ACV'*bCV - alphaMat(:,nodeIndex) + rho*(numNeighbor*xHatMatrixOld(:,nodeIndex) + sum(neighbor.*xHatMatrixOld,2)));
        end
        xHatMatrixNew = xHatMatrixNew2;
    else
        for ii = 1:numNode
            xHatMatrixOld = xHatMatrixNew;
            nodeInd = ii;
            nodeIndex = nodeInd;
            numNeighbor = sum(adjMatrixOrig(nodeInd,:));
            lastRCVDMat(:,:,nodeInd) = repmat(adjMatrixOrig(nodeInd,:),size(xHatMatrixNew,1),1).*xHatMatrixOld;
            alphaMat(:,nodeInd) = alphaMat(:,nodeInd) + rho*(numNeighbor*xHatMatrixOld(:,nodeInd) - sum(lastRCVDMat(:,:,nodeInd),2)) ;
            ACV = A(:,:,nodeIndex);
            bCV = b(:,nodeIndex);
            xHatMatrixNew2(:,nodeIndex) = inv(ACV'*ACV + 2*rho*numNeighbor*eye(numberVariable))*( ACV'*bCV - alphaMat(:,nodeIndex) + rho*(numNeighbor*xHatMatrixOld(:,nodeIndex) + sum(lastRCVDMat(:,:,nodeIndex),2)));
        end
        xHatMatrixNew = xHatMatrixNew2;
    end
    errorNormalized(jj,:) = sqrt(sum(sum((xHatMatrixNew - repmat(xHatCen,1,numNode)).^2)))/sqrt(sum(sum(repmat(xHatCen,1,numNode).^2)));  
    
end