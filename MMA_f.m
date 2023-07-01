function [objK,derivRadK] = MMA_f(x,L,U,p,q,r,j)
objK = sum((p(:,j)./(U-x)+q(:,j)./(x-L))) + r(j);
derivRadK = (p(:,j).*(U-x).^-2) - (q(:,j).*(x-L).^-2);
if 0
    numDeriv = zeros(size(derivRadK));
    eps = 0.001;
    for i = 1:length(x)
        temp = x;
        temp(i) = temp(i)+eps;
        numObj = sum((p(:,j)./(U-temp)+q(:,j)./(temp-L))) + r(j);
        numDeriv(i) = (numObj-objK)/eps;
    end
end
end
