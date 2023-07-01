function P = calcForces(nodeCord,forceBCs)

P = zeros(length(nodeCord)*2,1);

[len,~] = size(forceBCs);

for i =1:len %iterate on each force
    index1 = forceBCs(i,1)*2-1; %x direction
    index2 = forceBCs(i,1)*2; %y direction
    if forceBCs(i,3) == 1 %Outward pressure
        v = nodeCord(forceBCs(i,1),:);
        v = (v/norm(v))*forceBCs(i,2);
        P(index1) = v(1);
        P(index2) = v(2);
        %Need to adjust for node placement
    elseif forceBCs(i,3) == 0
        %No force
        5;
    elseif forceBCs(i,3) == 2
        v = [nodeCord(forceBCs(i,1),:) 0];
        v = cross(v,[0,0,1]);
        v = v(1:2);
        v = (v/norm(v))*forceBCs(i,2);
        P(index1) = v(1);
        P(index2) = v(2);
        
    elseif forceBCs(i,3) == 5
        P(index1) = forceBCs(i,2);
        P(index2) = 0;
        %Need to adjust for node placement
    elseif forceBCs(i,3) == 6
        P(index1) = 0;
        P(index2) = forceBCs(i,2);
        %Need to adjust for node placement
    end
    
end