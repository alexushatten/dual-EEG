function Y = scale_mat_0_1( X )
% scale_0_1: scale (aka normalize) values in matrix such that minimum = 0
% and maximum = 1
%
% Y = scale_0_1( X )
%
%  Inputs
% --------
% X: [n x m] numeric array of values
%
%  Outputs
% ---------
% Y: [n x m] numeric array of values, such that Y = X - min(X)/(max(X) - min(X)),
%            so that min(Y) = 0 and max(Y) = 1
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, May 2019
%-------------------------------------------------------------------------

Y = (X - min(X(:)))/(max(X(:))-min(X(:)));

end

