function X2Y = scale_X2Y( X, Y, ArrayWise )
% scale_X2Y: scale array X such that it has the same minimum and maximum as
% array Y
%
% X2Y = scale_X2Y( X, Y )
%
%  Inputs
% --------
% X: array
% Y: array
% ArrayWise: logical, whether to apply scaling to each column of X if X is
%            a matrix (= false, default), or whether to apply scaling across the
%            whole array instead (= true)
%
%  Outputs
% ---------
% X2Y: array X whose minimum and maximum equal minimum and maximum of Y
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, June 2019
%-------------------------------------------------------------------------
if nargin < 3
    ArrayWise = false;
end

Xs = scale_0_1(X, ArrayWise);
MinY = min(Y(:));
MaxY = max(Y(:));

X2Y = Xs*(MaxY-MinY)+MinY;

end

