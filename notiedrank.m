function R = notiedrank( X )
% NOTIEDRANK: compute the ranks of a sample, without adjusting for ties.
% Actually, in contrast to ranks without adjustment for ties, and
% contrarily to TIEDRANK, this function will assign each unique value of
% input X a non-unique integer value, like if you were converting the
% values to monotic integer values with equal increments, such that, for
% example:
%       [ 0.5, 0.5, 0.5,   1,   1,  10,  10,  10,  10 ]
% will become:
%       [   1,   1,   1,   2,   2,   3,   3,   3,   3 ]
%
% Therefore, this function does some kind of categorisation of values.
%
% R = notiedrank( X )
%
%  Inputs
% --------
% X: vector or matrix of values
%
%  Outputs
% ---------
% R: ranks of X
%
% See also TIEDRANK
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab&LREN, October 2019
%-------------------------------------------------------------------------

% [~,r] = sort(X);
% R = nan(size(X));
% R(r) = 1:numel(X);

U = unique(X);
R = nan(size(X));
Count = 0;
for u = 1:length(U)
    Count = Count+1;
    R(X==U(u))=Count;
end

% [~,I] = sort(X);
% [~,R] = sort(I);

end

