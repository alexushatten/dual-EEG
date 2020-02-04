function [Idx,Val,Median] = find_median( X )
% find_median: find value in vector closest to median
%
% [Idx,Val,Median] = find_median( X )
%
%  Inputs
% --------
% X: vector of values
%
%  Outputs
% ---------
% Idx: index of value closest to median in vector X
% Val: value in vector X closest to median
% Median: median value in vector X
% 
%
%-------------------------------------------------------------------------
% Cartool: https://sites.google.com/site/cartoolcommunity
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, November 2018
%-------------------------------------------------------------------------

Median = median(X);
[Val,Idx] = min(abs(X-Median));

end

