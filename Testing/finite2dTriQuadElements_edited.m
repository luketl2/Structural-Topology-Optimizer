%% Perform a 2D FEA of quadratic triangular element domain.
%Sunil Bean, University of Illinois, Argonne National Lab, 8/6/2020
%-------------------------------------------------------------------------
%This code assumes an isotropic material
%E:             Youngs Modulus, MPa
%v:             Poisson's Ratio
%t:             Material thickness, units mm
%nodeCord:      The coordinates of every node in the FEA domain, mm
%elemNodes:     For every element, the nodes of the element ordered 1,2,3,4,5,6
%constraint:    Nodes that are fixed BC (2 DOF fixed)
%BCs:           Any other boundary conditions or forces applied to the FEA
%-------------------------------------------------------------------------
%Outputs
%finalNodeCord: Coordinates after the FEA.Converted to mm.
%vmStress:      Von mises stress for each element. Converted to MPa
%dispU:         Nodal displacements from the FEA. Converted to mm.
%P:             Forces acting within the domain
%SM:            This is the inner stress matrix needed to calculate
                %stresses. It may make sense in future code to pass this as
                %a global variable to save computation time.
%Global Variables
%V:             Cholesky factorization of the free free stifness matrix
%Kff:           The free free part of the global stiffness matrix
%Kpf:           The forced free part of the global stiffness matrix
%fDOF:          The free degrees of freedom
%pDOF:          The prescribed degrees of freedom
%cDOF:          The user defined constrainted degrees of freedom
%KLocalDeriv:   The derivative of the local stiffness matrix wrt roe
%pen:           The penalty method exponent factor
function [finalNodeCord,vmStress,dispU,P,SM] = finite2dTriQuadElements_edited(E,v,t,nodeCord,elemNodes,elemArea,constraint,prescBCs,Pfree,roe)
overall_start = tic
global V
global Kff
global Kpf
global fDOF
global pDOF
global cDOF
global kLocalDeriv
global pen

pen = 4; % added here since the app cannot contain global vars
E = E*1e9; %Convert to Pa from GPa
t = t/1000; %Convert to m
nodeCord = nodeCord/1000; %Convert to m
%converting prescribed displacements to m
if ~isempty(prescBCs)
    [~,col] = size(prescBCs);
    if col ==1
        5;
    end
    for i = 2:col
        prescBCs(:,i) = prescBCs(:,i)/1000;
    end
end
%determine the number of elements
[elemN,~] = size(elemNodes);
%Initialze the final DOF displacement vector
dispU = zeros(length(nodeCord)*2,1);
%Initialize the degrees of freedom
dof = 2*length(nodeCord);
%Initialize the constitutive relation matrix, plane stress
E = (E/(1-v^2))*    [   1 v 0; ...
                        v 1 0; ...
                        0 0 (1-v)/2];
%Initialize the von mises constant matrix
VM = [      1        -1/2       0;...
            -1/2      1         0;...
            0         0         3];
%% Build the global stiffness matrix, initialize other matrices
%Initialize the global stiffness matrix
globalK = zeros(dof,dof);
%% Build local stiffness elements and assemble into global stiffness matrix
%First, there are shape function values that can be pre-calculated as they
%are the same for every element. We are applying the same Gauss Quadrature
%integration on every element. The partial derivativies of the shape
%functions with respect to their parent element coordinates is always the
%same at the gauss quadrature points. There will be 3 sets of partial
%dervatives to use.
weight = 1/6;
dNi1 = zeros(2,6);
dNi2 = zeros(2,6);
dNi3 = zeros(2,6);
% dNi4 = zeros(2,6); (removed)

for i = 1:6
    [dShape,~] = TriShape6Func(i,1/2,1/2);
    dNi1(:,i) = dShape;
    [dShape,~] = TriShape6Func(i,0,1/2);
    dNi2(:,i) = dShape;
    [dShape,~] = TriShape6Func(i,1/2,0);
    dNi3(:,i) = dShape;
    %     [dShape,~] = TriShape6Func(i,1/3,1/3);
    %     dNi4(:,i) = dShape; %triangle centroid for stress calculation
    %     (removed)
end

loop_1_start = tic
kLocalDeriv = zeros([12,12,elemN]);
for i = 1:elemN
    %Call upon the script to calculate local B matrix information. There
    %is 1 B matrix for each qudrature point.
    localIndex = elemNodes(i,:);
    localCoord = nodeCord(localIndex,:);
    [B1,B2,B3,detJ1,detJ2,detJ3] = localBtriQuad(dNi1,dNi2,dNi3,localCoord);
    %Now that the local B matrix is built at gauss points, calculate local
    %stiffness matrix. Density and penalty apply here.
    kLocal = weight*t*(B1'*E*(roe(i)^pen)*B1*detJ1+...
        B2'*E*(roe(i)^pen)*B2*detJ2+...
        B3'*E*(roe(i)^pen)*B3*detJ3);
    
    %Assemble the local into the global matrix
    localIndex = [(2*localIndex-1);(2*localIndex)];
    localIndex = localIndex(:);
    globalK(localIndex,localIndex) = globalK(localIndex,localIndex) + kLocal;
    
    %Now we want to calculate the derivative of local stiffness matrices
    %with respect to roe. This is used for sensitivities. 
    if pen > 0
        derivativeRoe = pen*roe(i)^(pen-1);
%     else
%         derivativeRoe = (1/(1-pen*(1-roe(i))))+(-pen)*(roe(i)/(1-pen*(1-roe(i)))^2);
    end
    kLocalDeriv(:,:,i) = weight*t* (B1'*E*derivativeRoe*B1*detJ1+...
        B2'*E*derivativeRoe*B2*detJ2+...
        B3'*E*derivativeRoe*B3*detJ3);
end
loop_1_time = toc(loop_1_start);

%% Apply BC's and solve for nodal displacements
%Eliminate contrained or prescribed DoF
cDOF = sort([constraint*2-1;constraint*2]);
if ~isempty(prescBCs) % if there are prescribed BC displacements
    uDOF = sort([prescBCs(:,1)*2-1;prescBCs(:,1)*2]);
else
    uDOF = [];
end
[pDOF,I] = sort([cDOF;uDOF]);
fDOF = [1:1:length(globalK)]'; % making a vector of the right length
fDOF(pDOF) = []; % removing the indices that are in pDOF
Pf = Pfree;
Pf(pDOF) = []; %Eliminate the prescribed indices
%need to combined constrained and prescribed into one Dp vector
if ~isempty(uDOF)
    Dp = [(prescBCs(:,2))';(prescBCs(:,3))'];
    Dp = Dp(:); %This is the prescribed component
else
    Dp = [];
end
Dp = [zeros(length(cDOF),1);Dp];
Dp = Dp(I);

%% Slove for nodal free free displacements, assemble into global 
% displacement vector
Kff = globalK(fDOF,fDOF);
Kfp = globalK(fDOF,pDOF);
Kpf = globalK(pDOF,fDOF);
Kpp = globalK(pDOF,pDOF);
V = chol(Kff,'lower');
Df = V'\(V\(Pf-Kfp*Dp)); % V*V'*Df = Pf-Kfp*Dp
nodeNum = 1:1:length(dispU);
dispU(setdiff(nodeNum,pDOF)) = Df;
dispU(pDOF) = Dp;

%% Solve for the forces
P = globalK*dispU;
%% Add the displacements to the current nodal coordinates
Dx = dispU(1:2:end);
Dy = dispU(2:2:end);
finalNodeCord = nodeCord+[Dx Dy];
%% Calculate the von Mises stresses, Output in MPa
vmStress = zeros(elemN,1);
SM = zeros(12,12,elemN);

loop_2_start = tic
for i = 1:elemN
    localIndex = [elemNodes(i,:)*2-1;elemNodes(i,:)*2];
    localIndex = localIndex(:);
    localIndex2 = elemNodes(i,:);
    localCoord = nodeCord(localIndex2,:);
    [B1,B2,B3,detJ1,detJ2,detJ3] = localBtriQuad(dNi1,dNi2,dNi3,localCoord);
    S = weight*(B1'*E*VM*E*B1*detJ1+B2'*E*VM*E*B2*detJ2+B3'*E*VM*E*B3*detJ3)/elemArea(i);
    SM(:,:,i) = S;
    vmStress(i) = sqrt((dispU(localIndex)'*S*dispU(localIndex)))/1e6;
end
loop_2_time = toc(loop_2_start);

%% Final conversions for export
finalNodeCord = finalNodeCord*1000;     %convert to mm
dispU = dispU*1000;    %convert to mm
end
