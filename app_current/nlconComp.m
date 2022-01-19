function [c,ceq,gc,gceq] = nlconComp(E,v,t,initialNodeCord,elemNodes,elemArea,constraint,prescBCs,Pfree,roe,objNodesIndices,objDir,initEnergy)
global Kff
global fDOF
global pDOF
global kLocalDeriv

ceq = 0; % always satisfied bc there is no eq constr
gceq = sparse(size(roe)); % apparently more efficient than zeros

%% adjust Pfree to have forces in the other direction to run FEA for constraint -- XXXXX maybe not right XXXXX
if mod(find(Pfree,1),2) % will be true if load is along x-axis, else it'll be along y
    Pfree = [0;Pfree]; % shift all the indices up by 1 to simulate y-axis force
    Pfree(end) = [];
else % adjust Pfree so FEA runs for a force along x-axis instead of y
    Pfree(1) = []; % shift all the indices down by 1
    Pfree(end+1) = 0; % add a 0 at the final y DOF
end
%% Performing FEA
[~,~,DispU,P,~] = finite2dTriQuadElements(E,v,t,initialNodeCord,elemNodes,elemArea,constraint,prescBCs,Pfree,roe);
DispU = DispU/1000; %convert to meters

%% calc inequality
if objDir == "y-axis" % this determines stiffness direction
    indices = zeros(size(objNodesIndices));
    for i = 1:length(objNodesIndices)
        indices(i) = objNodesIndices(i)*2-1;
    end
elseif objDir == "x-axis"
    indices = zeros(size(objNodesIndices));
    for i = 1:length(objNodesIndices)
        indices(i) = objNodesIndices(i)*2;
    end
else
    indices = zeros(2*length(objNodesIndices),1);
    for i = 1:length(objNodesIndices)
        indices(2*i-1) = objNodesIndices(i)*2-1;
        indices(2*i) = objNodesIndices(i)*2;
    end
end
c = P(indices)'*DispU(indices)-initEnergy*0.1; % constant sets max allowed energy
                                    % of objNode disp with force applied in N*m

%% calc constr function deriv wrt roh
[elemN,~] = size(elemNodes);
if nargout >1 % if the deriv is asked for
    PsiP = zeros(length(pDOF),1);
    PsiF = (-Kff^-1)'*P(fDOF);
    AdjointVector = [PsiP;PsiF];
    [~,I] = sort([pDOF;fDOF]);
    AdjointVector = AdjointVector(I);
    
    gc = zeros(elemN,1);
    for i = 1:elemN
        localIndex = [elemNodes(i,:)*2-1;elemNodes(i,:)*2];
        localIndex = localIndex(:);
        gc(i) = (AdjointVector(localIndex)'*kLocalDeriv(:,:,i)*DispU(localIndex));
    end
end