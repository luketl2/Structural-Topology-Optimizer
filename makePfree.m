% Makes an n*2 x 1 list of the load on each DOF of each node
% It is nx2 long bc each node has an x followed by a y DOF

function P = makePfree(nodeCoords,forceBCs,forceDir)

value = forceDir;
P = zeros(length(nodeCoords)*2,1); % preallocate size
if (value == "Up")
    for i = 1:length(forceBCs(:,1))
        P(forceBCs(i,1)*2) = forceBCs(i,2);
    end
elseif (value == "Down")
    for i = 1:length(forceBCs(:,1))
        P(forceBCs(i,1)*2) = -forceBCs(i,2);
    end
elseif (value == "Left")
    for i = 1:length(forceBCs(:,1))
        P(forceBCs(i,1)*2-1) = -forceBCs(i,2);
    end
elseif (value == "Right")
    for i = 1:length(forceBCs(:,1))
        P(forceBCs(i,1)*2-1) = forceBCs(i,2);
    end
end

end

