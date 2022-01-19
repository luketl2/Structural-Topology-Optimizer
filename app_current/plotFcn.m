function stop = outFcn(x,optimValues,state)
switch state
    case 'init'    
        hold on  
        
    case 'iter'       
        % Concatenate current point and objective function
        % value with history. x must be a row vector.        
        history.fval = [history.fval; optimValues.fval];      
        history.x = [history.x; x];       
        % Concatenate current search direction with        
        % searchdir.        
        delete(findobj('type', 'patch'));        
        for i = 1:elemN           
            nodes = elemNodes(:,1:3);           
            xy = initialNodeCord(nodes(i,1:3),:);           
            patch(xy(:,1),xy(:,2),x(i))           
        end       
        pause(.01)        
        colorbar        
        caxis([0 1])        
        axis equal        
        axis off
        view(2)

    case 'done'
        hold off
        
    otherwise    
end
stop = false;
end
