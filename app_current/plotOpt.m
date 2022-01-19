function plotOpt(elemNodes,finalNodeCord,roe,ax)

%determine the number of elements
[elemN,~] = size(elemNodes);

%% Plotting of Results
for i = 1:elemN
    nodes = elemNodes(:,1:3);
    xy = finalNodeCord(nodes(i,1:3),:);
    patch(ax,xy(:,1),xy(:,2),roe(i))
end
colorbar(ax)
axes(ax)
axis equal
view(2)
end