function plotFEA(vmStress,elemNodes,finalNodeCord,roe,ax,maxVM)
if max(roe) > 0.95
    maxVM.Value = (max(vmStress(roe>0.95)));
    display(maxVM.Value)
    %determine the number of elements
    [elemN,~] = size(vmStress);
    caxis(ax,[0,maxVM.Value])

    %% Plotting of Results
    vmStressNorm = vmStress;
    for i = 1:elemN
        nodes = elemNodes(:,1:3);
        xy = finalNodeCord(nodes(i,1:3),:);
        patch(ax,xy(:,1),xy(:,2),vmStressNorm(i))
    end
    colormap(ax,jet);
    colorbar(ax)
    axes(ax)
    axis equal
    view(2)
end
end