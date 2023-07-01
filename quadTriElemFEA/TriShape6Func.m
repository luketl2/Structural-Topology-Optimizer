%% Output of a 6 node triangular isoparametric element
function [dShape,Shape] = TriShape6Func(FunNum,Z,E)
%Inpute the desired shape function and derivatives
%The partial derivative of Z comes first out of the dShape method
if FunNum == 1
    dShape = [4*Z-1;0]; %fixed
    Shape = 2*Z*(Z-1/2); 
elseif FunNum == 2
    dShape = [0;4*E-1]; %fixed
    Shape = 2*E*(E-1/2); 
elseif FunNum == 3
    dShape = [4*Z+4*E-3;4*E+4*Z-3]; %fixed
    Shape = 2*(1-E-Z)*(0.5-E-Z); %this isnt right, 0.5 @ 0,0
elseif FunNum == 4
    dShape = [4*E;4*Z]; %fixed
    Shape = 4*E*Z;
elseif FunNum == 5
    dShape = [-4*E;4-8*E-4*Z]; %fixed
    Shape = 4*E*(1-E-Z);
elseif FunNum == 6
    dShape = [4-4*E-8*Z;-4*Z]; %fixed
    Shape = 4*Z*(1-E-Z);
else
    error('Input a shape function number from 1 to 6');
end

end