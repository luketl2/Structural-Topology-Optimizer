function [objK,derivRadK] = ComplianceSensitivity(E,v,t,initialNodeCord,elemNodes,elemArea,constraint,prescBCs,Pfree,roe,objNodesIndices,objDir)
global Kff
global fDOF
global pDOF
global kLocalDeriv

%% Performing FEA
[~,~,DispU,P,~] = finite2dTriQuadElements(E,v,t,initialNodeCord,elemNodes,elemArea,constraint,prescBCs,Pfree,roe);
DispU = DispU/1000; %convert to meters

%% calc obj function
if objDir == "x-axis" % this determines stiffness direction
    indices = 2*objNodesIndices - 1;
elseif objDir == "y-axis"
    indices = 2*objNodesIndices;
else % objDir == "both"
    indices = [2*objNodesIndices - 1; 2*objNodesIndices]; % order doesn't matter for dot product
end
objK = -(P(indices)'*DispU(indices)); % so we dot the force in the DOFs of each objNode w/ the disp of each DOF for all nodes

%% calc obj function deriv wrt roh -- XXXXXX might need to change this section depending on obj function
[elemN,~] = size(elemNodes);
if nargout >1 % if the deriv is asked for
    PsiP = zeros(length(pDOF),1);
    PsiF = (-Kff^-1)'*P(fDOF);
    AdjointVector = [PsiP;PsiF];
    [~,I] = sort([pDOF;fDOF]);
    AdjointVector = AdjointVector(I);
    
    derivRadK = zeros(elemN,1);
    for i = 1:elemN
        localIndex = [elemNodes(i,:)*2-1;elemNodes(i,:)*2];
        localIndex = localIndex(:);
        derivRadK(i) = -(AdjointVector(localIndex)'*kLocalDeriv(:,:,i)*DispU(localIndex));
    end
end

