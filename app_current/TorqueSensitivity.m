function [objK,derivRadK] = TorqueSensitivity(E,v,t,initialNodeCord,elemNodes,elemArea,constraint,prescBCs,Pfree,roe,desiredK,MaxTorqueK,angRange,filtW)
global V
global Kff %#ok<NUSED>
global Kpf
global fDOF
global pDOF
global kLocalDeriv
global Lrot
global pen

%% Calculating derivative of roe
[elemN,~] = size(elemNodes);
%% Performing FEA
plotTrue = false;
[~,~,DispU,P,~] = finite2dTriQuadElements(E,v,t,initialNodeCord,elemNodes,elemArea,constraint,prescBCs,Pfree,filtW*roe,plotTrue);
DispU = DispU/1000; %convert to meters

%% Calculate the radial stiffness
% objK = abs((L'*P)/(angRange*(pi/180)))/MaxTorqueK - desiredK;
objK = ((Lrot'*P) - desiredK)/MaxTorqueK;
% objK = abs((L'*P)/(angRange*(pi/180)))/abs(((L'*Pang)/(angRange*(pi/180)))) - desiredK;

%% Calcuate the radial stifness wrt derivative of Roe 
if nargout >1 % if the deriv is asked for
    PsiP = Lrot(pDOF);
    PsiF = -V'\(V\(Lrot(pDOF)'*Kpf)');
%     PsiF = -V\(V'\(Kpf'*Lrot(pDOF)));
    AdjointVector = [PsiP;PsiF];
    [~,I] = sort([pDOF;fDOF]);
    AdjointVector = AdjointVector(I);
    
    derivRadK = elemN;
    for i = 1:elemN
        localIndex = [elemNodes(i,:)*2-1;elemNodes(i,:)*2];
        localIndex = localIndex(:);
        derivRadK(i) = (AdjointVector(localIndex)'*kLocalDeriv(:,:,i)*DispU(localIndex))/MaxTorqueK;
    end
    
    derivRadK = (derivRadK*filtW')';
end

end