function [Order, R] = angular_order_eig( R )
% angular_order_eig: find order (and re-order) of variables in correlation
% matrix based on AOE (angular order of eigenvectors)
%
% [Order, R] = angular_order_eig( R )
%
%  Inputs
% --------
% R: correlation matrix
%
%  Outputs
% ---------
% Order: order of columns based on angular order of eigenvectors
%
%-------------------------------------------------------------------------
% NB: based on R package "corrplot"
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, November 2018
%-------------------------------------------------------------------------

[n,m] = size(R);
if n ~= m
    error('R should be square')
end
[V,~] = eig(R);
e1 = V(:,1);
e2 = V(:,2);
if all(e1>0)
    Alpha = atan(e2 ./ e1);
else
    Alpha = atan(e2 ./ e1) + pi;
end
[~,Order] = sort(Alpha);
if nargout>1
    R = R(Order,Order);
end

end

