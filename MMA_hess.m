function h = MMA_hess(x,lambda,L,U,p,q)
h = diag(2.*(p(:,1).*(U-x).^-3) + 2.*(q(:,1).*(x-L).^-3));
if length(p(1,:)) > 1
    for j = 2:length(p(1,:))
        h = h + lambda.ineqnonlin(j-1)*diag(2.*(p(:,1).*(U-x).^-3) + 2.*(q(:,1).*(x-L).^-3));
    end
end

