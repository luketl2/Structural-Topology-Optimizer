function [objK,derivRadK] = StiffnessSensitivity(E,v,t,initialNodeCord,elemNodes,elemArea,constraint,prescBCs,Pfree,roe,objNodesIndices,objDir)
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

objK = P(indices)'*DispU(indices); % only calcs for nodes objNodes

%% calc obj function deriv wrt roe
[elemN,~] = size(elemNodes);
if nargout > 1 % if the deriv is asked for
    PsiP = zeros(length(pDOF),1);
    PsiF = (-Kff^-1)'*P(fDOF);
    AdjointVector = [PsiP;PsiF];
    [~,I] = sort([pDOF;fDOF]);
    AdjointVector = AdjointVector(I);
    
    derivRadK = zeros(elemN,1);
    for i = 1:elemN
        localIndex = [elemNodes(i,:)*2-1;elemNodes(i,:)*2];
        localIndex = localIndex(:);
        derivRadK(i) = (AdjointVector(localIndex)'*kLocalDeriv(:,:,i)*DispU(localIndex));
    end
    
% Tried without loop but this einsum is different than numpy.einsum
%     addpath(strcat(cd,'\einsum'))
%     temp = elemNodes';
%     temp = [temp(:)*2 - 1, temp(:)*2]';
%     temp = temp(:);
%     local_indices = reshape(temp, [12, elemN]);
%     temp = repmat(AdjointVector, [1, elemN]);
%     adjoint_matrix = temp(local_indices);
%     temp = repmat(DispU, [1, elemN]);
%     disp_matrix = temp(local_indices);
%     p1 = einsum(adjoint_matrix, kLocalDeriv, 'il, ijk -> jk')
%     derivRadK = einsum(p1, disp_matrix, 'il, ik -> k')
    
end

