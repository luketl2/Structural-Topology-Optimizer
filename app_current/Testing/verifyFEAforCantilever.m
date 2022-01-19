%This code is to test the accuracy of the stress calculation for the FEA
%program
close all
%The path is added for the location of the distmesh, which is the
%Persson Per-Olof & Strang, Gilbert code for mesh generation.
addpath(strcat(cd,'\distmesh'))
%The path for the FEA program
addpath(strcat(cd,'\quadTriElemFEA'))
%% [1] Initial Conditions for the particular FEA problem - ASTM A-36 steel
global pen 
E =         211;    %Youngs Modulus, GPa
v =         0.29;   %Poissons ration, unitless
t =         10;     %Material thickness, mm
scale =     100;    %mm (how long is half a side of the plate)
BCtype =    6;      %Boundary condition type, 0 is no presribed force, 
                    %1 is outward pressure, 2 is radial clockwise, 5 right 
                    %pressure, 6 upward pressure (to be made more elegant
                    %later)
pen =       3;      %penalty method exponent, null
FORCE =     1000;   %Newtons

%% [2] Creating a mesh, this is based on Persson Per-Olof & Strang, Gilbert 
%code for mesh generation Example: (Rectangle)
fd=@(p) drectangle(p,-1,1,-.5,.5);
pfix=[-1,-.5;-1,.5;1,-.5;1,.5];
[nodeCord,elemNodes]=distmesh2d(fd,@huniform,0.1,[-1,-.5;1,.5],pfix);

%% [3] Adding FEA mid nodes, the Strang code does not do this.
[elemNodes,nodeCord] = addMidNodes(nodeCord,elemNodes);

%% [4] Initialize the density matrix based on the amount of elements, give 
%every element a density of 1
[elemN,~] = size(elemNodes);
roe = ones(elemN,1);

%% [5] Find the edge nodes for constraints and force application 
%respectively this is problem specific, because the mesh is initially 
%normilized to a length scale of +/- 1, it is easy to identify the edges.
bottomNodes = find(nodeCord(:,2)<= -.4999);
topNodes = find(nodeCord(:,2)>= .4999);
leftNodes = find(nodeCord(:,1)<= -.9999);
rightNodes = find(nodeCord(:,1)>= .9999);

%% [6] Define constraints
constraint = leftNodes;%These are nodes that are fixed.
uPrescribed = []; %There are no prescribed displacements for this problem

%% [7] Initialize the Boundary Conditions
% this code is for symmetric force on right side
% forceDir = "Down";
% forceSize = ceil(length(rightNodes)/10);
% forceBCs = zeros(forceSize,2);
% forceBCs = forceBCs + FORCE/forceSize;
% rightNodes(ceil(length(rightNodes)/2-(forceSize)/2):floor(length(rightNodes)/2+(forceSize)/2));

% this is for point force at top right corner
% [forceNodes,~] = intersect(rightNodes,topNodes);
% forceBCs = zeros(length(forceNodes)) + FORCE/length(forceNodes);

% this code is for force along entire right side of beam
forceDir = "Down";
forceNodes = rightNodes;
forceBCs = zeros(size(forceNodes)) + FORCE/length(forceNodes);

forceBCs = [forceNodes,forceBCs,zeros(size(forceNodes))+BCtype]; % [node number, total force on the node, type of force (helps determine force components)]
Pfree = makePfree(nodeCord,forceBCs,forceDir);
prescBCs = [];

%% [8] Scale the nodal coordinates
nodeCord = (nodeCord* scale);

%% [9] Calculate the area of each element, m^2
elemArea = zeros(length(elemNodes),1);
for i = 1:length(elemNodes)
    n1 = nodeCord(elemNodes(i,1),:);
    n2 = nodeCord(elemNodes(i,2),:);
    n3 = nodeCord(elemNodes(i,3),:);
    x = [n1(1),n2(1),n3(1)]/1000; %convert x cord to m
    y = [n1(2),n2(2),n3(2)]/1000; %convert y cord to m
    elemArea(i) = polyarea(x,y); %here the integral is over the area, so 
                                %"volume" is really just area for 2D
end

% %% [10] Set the desired plotting
% plotTrue = 1;
% 
%% [11] Performing the FEA
[finalNodeCord,vmStress,disp,P] = finite2dTriQuadElements(E,v,t,nodeCord,elemNodes,elemArea,constraint,prescBCs,Pfree,roe);
disp = disp/1000;

%% [12] Using this for checking the sensitivity function objDeriv calcs
objNodeIndices = rightNodes; % for the cantilever example
numericalObjDeriv = numericalDerivCalc(E,v,t,nodeCord,elemNodes,elemArea,constraint,prescBCs,Pfree,roe,objNodeIndices);
[~,adjointObjDeriv] = StiffnessSensitivity(E,v,t,nodeCord,elemNodes,elemArea,constraint,prescBCs,Pfree,roe,objNodeIndices,"y-axis");
error = numericalObjDeriv-adjointObjDeriv;
totalError = sum(abs(error))

%% [13] Plotting
figure('DefaultAxesFontSize',14,'DefaultAxesFontName',"Times New Roman")
hold on
title("Objective Sensitivities for Numerical and Adjoint Methods","FontSize",20)
xlabel("element number")
ylabel("objective sensitivity (J)")
plot(adjointObjDeriv,'--',"DisplayName","Adjoint Method", "LineWidth", 2)
plot(numericalObjDeriv1,'-.',"DisplayName","Numerical Method, Step-size: 0.1", "LineWidth", 2)
plot(numericalObjDeriv5,'-s',"DisplayName","Numerical Method, Step-size: 0.05", "LineWidth", 2)
plot(numericalObjDeriv2,'-*',"DisplayName","Numerical Method, Step-size: 0.02", "LineWidth", 2)
legend('Location','se')
xlim([200,210]);
hold off
