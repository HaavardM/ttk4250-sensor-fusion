function jacMatrix = jacobianFD(fun,x,dx)

% General finite-difference Jacobian
% Both fun and x are assumed vector valued

f = fun(x);
m = size(x,1);
n = size(f,1);

jacMatrix = zeros(n,m);
for k=1:size(x,1)
    xPerturbed = x;
    xPerturbed(k) = xPerturbed(k) + dx;
    fPerturbed = fun(xPerturbed);
    jacMatrix(:,k) = (fPerturbed-f)/dx;
end
