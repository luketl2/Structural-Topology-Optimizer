function [B1,B2,B3,detJ1,detJ2,detJ3] = localBtriQuad(dNi1,dNi2,dNi3,localCoord)
%Determine Jacobian transformations and B matricies at gauss quadrature
%ponts

Jac1 = (dNi1*localCoord);
Jac2 = (dNi2*localCoord);
Jac3 = (dNi3*localCoord);

P1 = Jac1\dNi1;
P2 = Jac2\dNi2;
P3 = Jac3\dNi3;

B1 = [P1(1,1),0,P1(1,2),0,P1(1,3),0,P1(1,4),0,P1(1,5),0,P1(1,6),0;...
      0,P1(2,1),0,P1(2,2),0,P1(2,3),0,P1(2,4),0,P1(2,5),0,P1(2,6);...
      P1(2,1),P1(1,1),P1(2,2),P1(1,2),P1(2,3),P1(1,3),P1(2,4),P1(1,4),P1(2,5),P1(1,5),P1(2,6),P1(1,6)];

B2 = [P2(1,1),0,P2(1,2),0,P2(1,3),0,P2(1,4),0,P2(1,5),0,P2(1,6),0;...
      0,P2(2,1),0,P2(2,2),0,P2(2,3),0,P2(2,4),0,P2(2,5),0,P2(2,6);...
      P2(2,1),P2(1,1),P2(2,2),P2(1,2),P2(2,3),P2(1,3),P2(2,4),P2(1,4),P2(2,5),P2(1,5),P2(2,6),P2(1,6)];
  
B3 = [P3(1,1),0,P3(1,2),0,P3(1,3),0,P3(1,4),0,P3(1,5),0,P3(1,6),0;...
      0,P3(2,1),0,P3(2,2),0,P3(2,3),0,P3(2,4),0,P3(2,5),0,P3(2,6);...
      P3(2,1),P3(1,1),P3(2,2),P3(1,2),P3(2,3),P3(1,3),P3(2,4),P3(1,4),P3(2,5),P3(1,5),P3(2,6),P3(1,6)];

detJ1 = det(Jac1);
detJ2 = det(Jac2);
detJ3 = det(Jac3);

end