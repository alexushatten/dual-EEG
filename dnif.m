function xdi = dnif(idx,L)
% DNIF (aka reverse find): get logical vector from subscript indices
% obtained with find.m
%
%   idx:    vector of subscript indices
%   L:      length of the logical vector requested
%
% Renaud Marquis @ FBM lab

xdi = false(L,1);
xdi(idx) = true;

end