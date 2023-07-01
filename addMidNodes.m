function [elemNodes,nodeCord] = addMidNodes(initCord,initNodes)

%find the index to start creating new nodes
newIndex = length(initCord)+1;
%extend coordinate matrix to accomodate half nodes
initNodes = [initNodes,zeros(size(initNodes))];

%iterate through each element and add new nodes
for i = 1:length(initNodes)
    nodes = initNodes(i,:);
    
    location1 = 4;
    %go from 1 --> 2
    if initNodes(i,location1) == 0 %Current element node needs to be filled in
        %Check for if there is another element side
        node1 = nodes(1);
        node2 = nodes(2);
        [x1,~] = find(initNodes == node1);
        [x2,~] = find(initNodes == node2);
        commonSideElements = intersect(x1,x2);
        
        if length(commonSideElements) > 1 %there is another elemnt w/common side
            %need to find the index in the other element this node should exist
            %as "location2". This assumes only one side is shared between two
            %elements
            
            element2 = commonSideElements(intersect(x1,x2)~= i);
            nodes2 = initNodes(element2,:);
            location2 = [find(nodes2 == node1),find(nodes2 == node2)];
            location2 = sort(location2);
            
            if all(location2 == [1,2])
                location2 = 4;
            elseif all(location2 == [2,3])
                location2 = 5;
            elseif all(location2 == [1,3])
                location2 = 6;
            end
            
            initNodes(i,location1) = newIndex;
            initNodes(element2,location2) = newIndex;
            
        else %This node is not shared but needs to be filled in
            initNodes(i,location1) = newIndex;
        end
        
        %Need to add the new node and location to the current nodal matrix
        initCord(newIndex,1) = (initCord(node1,1)+initCord(node2,1))/2;
        initCord(newIndex,2) = (initCord(node1,2)+initCord(node2,2))/2;
        
        newIndex = newIndex +1;
    end
    
    location1 = 5;
    %go from 2 --> 3
    if initNodes(i,location1) == 0 %Current element node needs to be filled in
        %Check for if there is another element side
        node1 = nodes(2);
        node2 = nodes(3);
        [x1,~] = find(initNodes == node1);
        [x2,~] = find(initNodes == node2);
        commonSideElements = intersect(x1,x2);
        
        if length(commonSideElements) > 1 %there is another elemnt w/common side
            %need to find the index in the other element this node should exist
            %as "location2". This assumes only one side is shared between two
            %elements
            
            element2 = commonSideElements(intersect(x1,x2)~= i);
            nodes2 = initNodes(element2,:);
            location2 = [find(nodes2 == node1),find(nodes2 == node2)];
            location2 = sort(location2);
            
            if all(location2 == [1,2])
                location2 = 4;
            elseif all(location2 == [2,3])
                location2 = 5;
            elseif all(location2 == [1,3])
                location2 = 6;
            end
            initNodes(i,location1) = newIndex;
            initNodes(element2,location2) = newIndex;
            
        else %This node is not shared but needs to be filled in
            initNodes(i,location1) = newIndex;
        end
        
        %Need to add the new node and location to the current nodal matrix
        initCord(newIndex,1) = (initCord(node1,1)+initCord(node2,1))/2;
        initCord(newIndex,2) = (initCord(node1,2)+initCord(node2,2))/2;
        
        newIndex = newIndex +1;
    end
    
    location1 = 6;
    %go from 3 --> 1
    if initNodes(i,location1) == 0 %Current element node needs to be filled in
        %Check for if there is another element side
        node1 = nodes(3);
        node2 = nodes(1);
        [x1,~] = find(initNodes == node1);
        [x2,~] = find(initNodes == node2);
        commonSideElements = intersect(x1,x2);
        
        if length(commonSideElements) > 1 %there is another elemnt w/common side
            %need to find the index in the other element this node should exist
            %as "location2". This assumes only one side is shared between two
            %elements
            
            element2 = commonSideElements(intersect(x1,x2)~= i);
            nodes2 = initNodes(element2,:);
            location2 = [find(nodes2 == node1),find(nodes2 == node2)];
            location2 = sort(location2);
            
            if all(location2 == [1,2])
                location2 = 4;
            elseif all(location2 == [2,3])
                location2 = 5;
            elseif all(location2 == [1,3])
                location2 = 6;
            end
            initNodes(i,location1) = newIndex;
            initNodes(element2,location2) = newIndex;
            
        else %This node is not shared but needs to be filled in
            initNodes(i,location1) = newIndex;
        end
        
        %Need to add the new node and location to the current nodal matrix
        initCord(newIndex,1) = (initCord(node1,1)+initCord(node2,1))/2;
        initCord(newIndex,2) = (initCord(node1,2)+initCord(node2,2))/2;
        
        newIndex = newIndex +1;
    end
end

elemNodes = initNodes;
nodeCord = initCord;

end