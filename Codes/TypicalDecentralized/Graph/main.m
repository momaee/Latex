%%
%%
clc
clear all
close all
tic
%% Parameters
%Network
numNode = 50;
compGraphNumEdge = numNode*(numNode-1)/2;
connectivityRatio = 0.5;
connectivityRatio = max([(numNode-1)/compGraphNumEdge, connectivityRatio]); 
numEdge = ceil(connectivityRatio*compGraphNumEdge);
treeFlag = 0;
branchFactor = 0;
if numEdge == numNode-1
    prompt = ['The Underlying Graph is a tree, Please enter the branch factor, it should be less than',num2str(numNode-1)];
    branchFactor = input(prompt);
    treeFlag = 1;
end
regularFlag = 0;
degree = 0;
if regularFlag == 1
    degree = ceil(2*numEdge/numNode);
    numEdge = degree*numNode/2;
end
connectivityRatio = numEdge/compGraphNumEdge;
%Local Objectives
numVariable = 3;
numObservation = 3;
noiseVariance = 1e-3;
%Simulation
numDiffNetwork = 1;
numExperiment = 1;
numIteration = 500;
errorNormalized = zeros(numIteration,numExperiment,numDiffNetwork);

%%
for netInedx = 1:numDiffNetwork
    %%%%%Network Generation
    [adjacencyMat] = networkGeneration(treeFlag,branchFactor,regularFlag,degree,numNode,numEdge);
    for expIndex = 1:numExperiment
        %%% Experiment Generation
        [x,observationVec,measurementMat,LipshitzCons,strongConv,AA,BB,mf,Mf] = experimetGeneraton(numNode,numVariable,numObservation,noiseVariance);
        %%% Ceneralized solution
        xHatCen = inv(AA'*AA)*AA'*BB;
        %%% ADMM Algorithm
        c = 0.1; %ADMMM Parameter        
        [errorNormalized(:,expIndex,netInedx)] = ADMM(adjacencyMat,c,numNode,numVariable,numIteration,measurementMat,observationVec,xHatCen);
           
    end
end
semilogy(mean(errorNormalized(:,:,1),2))
toc




