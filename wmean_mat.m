function WXY = wmean_mat( X, Y, W )
% wmean_mat: weighted average for matrices (and vectors)
%
% WXY = wmean_mat( X, Y, W )
%
%  Inputs
% --------
% X: [n x m] first matrix (can also be a vector)
% Y: [n x m] second matrix (can also be a vector)
% W: [1 x 2] vector with weights corresponding to X (W(1)) and Y (W(2))
%
%  Outputs
% ---------
% WXY: [n x m] weighted average matrix (or vector) between X and Y with
%      weights W
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, June 2019
%-------------------------------------------------------------------------

if ~all(size(X) == size(Y))
    error('Input matrices X and Y sizes do not match')
end

WXY = ( X * W(1) + Y * W(2) ) / sum(W);

end

