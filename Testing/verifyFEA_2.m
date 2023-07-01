%This code is to test the accuracy of the stress calculation for the FEA
%program
close all
%The path is added for the location of the distmesh, which is the
%Persson Per-Olof & Strang, Gilbert code for mesh generation.
addpath(strcat(cd,'\distmesh'))
%The path for the FEA program
addpath(strcat(cd,'\quadTriElemFEA'))
%% [1] Initial Conditions for the particular FEA problem
global pen 
E =         200;    %Youngs Modulus, GPa
v =         0.3;    %Poissons ration, unitless
t =         10;     %Material thickness, mm
scale =     100;    %mm (how long is half a side of the plate)
BCtype =    6;      %Boundary condition type, 0 is no presribed force, 
                    %1 is outward pressure, 2 is radial clockwise, 5 right 
                    %pressure, 6 upward pressure (to be made more elegant
                    %later)
pen =       3;      %penalty method exponent, null
FORCE =     1000;   %Newtons

%% [2] Creating a mesh, this is based on Persson Per-Olof & Strang, Gilbert 
%code for mesh generation Example: (Rectangle with circular hole, refined 
%at circle boundary)
fd=@(p) ddiff(drectangle(p,-1,1,-1,1),dcircle(p,0,0,0.5));
fh=@(p) 0.05+0.3*dcircle(p,0,0,0.02);
[nodeCord,elemNodes]=distmesh2d(fd,fh,0.04,[-1,-1;1,1],[-1,-1;-1,1;1,-1;1,1]);

%% [3] Adding FEA mid nodes, the Strang code does not do this.
[elemNodes,nodeCord] = addMidNodes(nodeCord,elemNodes);

%% [4] Initialize the density matrix based on the amount of elements, give 
%every element a density of 1
[elemN,~] = size(elemNodes);
roe = ones(elemN,1);

%% [5] Find the bottom & top nodes for constraints and force application 
%respectively this is problem specific, because the mesh is initially 
%normilized to a length scale of +/- 1, it is easy to identify the edges.
bottomNodes = find(nodeCord(:,2)<= -.98);
topNodes = find(nodeCord(:,2)>= .98);

%% [6] Define constraints
constraint = bottomNodes;%These are nodes that are fixed.
uPrescribed = []; %There are no prescribed displacements for this problem

%% [7] Initialize the Boundary Conditions
forceBCs = zeros(size(topNodes));
forceBCs = forceBCs + FORCE/length(topNodes); 
forceBCs = [topNodes,forceBCs,zeros(size(topNodes))+BCtype];
Pfree = calcForces(nodeCord,forceBCs);
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

%% [10] Set the desired plotting
plotTrue = 3;

%% [11] Performing the FEA
[finalNodeCord,~,disp,P] = finite2dTriQuadElements(E,v,t,nodeCord,elemNodes,elemArea,constraint,prescBCs,Pfree,roe,plotTrue);
disp = disp/1000;