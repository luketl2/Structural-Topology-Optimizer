function [numericalDeriv] = numericalDerivCalc(E,v,t,initialNodeCord,elemNodes,elemArea,constraint,prescBCs,Pfree,roe,objNodesIndices)
%% Performing initial FEA
[~,~,DispU,P,~] = finite2dTriQuadElements(E,v,t,initialNodeCord,elemNodes,elemArea,constraint,prescBCs,Pfree,roe);
DispU = DispU/1000; %convert to meters
    
%% calc initial obj function
indices = zeros(2*length(objNodesIndices),1);
for i = 1:length(objNodesIndices)
    indices(2*i) = objNodesIndices(i)*2;
    indices(2*i-1) = objNodesIndices(i)*2-1;
end
objK0 = P(indices)'*DispU(indices); % only calcs for nodes objNodes

%% perturb densities and re-run FEA and re-calc new objectives to find deriv vec
numericalDeriv = zeros(size(roe));
for i = 1:length(roe)
    temp = roe;
    temp(i) = temp(i)-0.02;
    [~,~,DispU,P,~] = finite2dTriQuadElements(E,v,t,initialNodeCord,elemNodes,elemArea,constraint,prescBCs,Pfree,temp);
    DispU = DispU/1000; %convert to meters
    
    objK = P(indices)'*DispU(indices); % only calcs for nodes objNodes
    numericalDeriv(i) = (objK-objK0)/(temp(i)-roe(i));
end
end
