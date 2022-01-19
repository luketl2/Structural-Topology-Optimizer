function plotFEA(plotTrue,vmStress,elemNodes,finalNodeCord,roe)

%determine the number of elements
[elemN,~] = size(vmStress);

%% Plotting of Results
if plotTrue == 1
    vmStressNorm = vmStress;
    figure
    hold on
    for i = 1:elemN
        nodes = elemNodes(:,1:3);
        xy = finalNodeCord(nodes(i,1:3),:);
        patch(xy(:,1),xy(:,2),vmStressNorm(i))
    end
    colorbar
    axis equal
    axis off
    view(2)
    hold off
elseif plotTrue == 2
    figure
    hold on
    for i = 1:elemN
        nodes = elemNodes(:,1:3);
        xy = finalNodeCord(nodes(i,1:3),:);
        patch(xy(:,1),xy(:,2),roe(i))
    end
    colorbar
    axis equal
    axis off
    view(2)
    hold off
elseif plotTrue == 3
    vmStressNorm = vmStress;       
    %     index  = (vmStressNorm >= 250);
    %     vmStressNorm(index) = 250;
    figure
    hold on
    for i = 1:elemN
        if roe(i) > 0.95
            nodes = elemNodes(:,1:3);
            xy = finalNodeCord(nodes(i,1:3),:);
            patch(xy(:,1),xy(:,2),vmStressNorm(i))
        end
    end
    colorbar
    axis equal
    axis off
    view(2)
    hold off
    figure
    hold on
    for i = 1:elemN
        nodes = elemNodes(:,1:3);
        xy = finalNodeCord(nodes(i,1:3),:);
        patch(xy(:,1),xy(:,2),roe(i))
    end
    colorbar
    axis equal
    axis off
    view(2)
    hold off
    
    %Output the max von mises on a more solid element
    display('The max von mises stress is:')
    display(max(vmStress(roe>0.95)))
end

end