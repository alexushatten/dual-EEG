function Y = scale_0_1( X, ArrayWise )
% scale_0_1: scale (aka normalize) values in vector such that minimum = 0
% and maximum = 1
%
% Y = scale_0_1( X )
%
%  Inputs
% --------
% X: [n x m] numeric array of values, scaling is applied for each column if X
%            is a matrix (m > 1)
% ArrayWise: logical, whether to apply scaling to each column of X if X is
%            a matrix (= false, default), or whether to apply scaling across the
%            whole array instead (= true)
%
%  Outputs
% ---------
% Y: [n x m] numeric array of values, such that Y = X - min(X)/(max(X) - min(X)),
%            so that min(Y) = 0 and max(Y) = 1
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, May 2019, last updated June 2019
%-------------------------------------------------------------------------

if nargin < 2
    ArrayWise = false;
end
if ~ArrayWise
    Min = min(X);
    Max = max(X);
else
    Min = min(X(:));
    Max = max(X(:));
end
if numel(X) == length(X) || ArrayWise
    if (Max-Min)==0
        warning('''X'' is constant and leads to division by 0, setting output to 0!')
        Y = zeros(size(X));
    else
        Y = (X - Min)/(Max-Min);
    end
else
    Y = bsxfun(@rdivide,bsxfun(@minus,X,Min),Max-Min);
    if any(Max==Min)
        warning('Some of the columns of ''X'' are constant and leads to division by 0, setting them to 0')
        Y(isnan(Y))=0;
    end
end

