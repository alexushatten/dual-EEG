function [Idx,Val,Mean] = find_mean( X )
% find_mean: find value in vector closest to mean
%
% [Idx,Val,Mean] = find_mean( X )
%
%  Inputs
% --------
% X: vector of values
%
%  Outputs
% ---------
% Idx: index of value closest to mean in vector X
% Val: value in vector X closest to mean
% Mean: mean value in vector X
% 
%
%-------------------------------------------------------------------------
% Cartool: https://sites.google.com/site/cartoolcommunity
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, November 2018
%-------------------------------------------------------------------------

Mean = mean(X);
[Val,Idx] = min(abs(X-Mean));

end

