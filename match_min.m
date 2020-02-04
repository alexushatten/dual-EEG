function [Idx, IdxVal] = match_min( X, Y, Function )
%MATCH_MIN: correspondance of X in Y based on smallest absolute value
%estimated based on function Function
%
% Idx = match_min( X, Y, Function )
%
%  Inputs
% --------
% X: [n x m] array of coordinates: row represent observations, columns
%            represent dimensions
% Y: [n x m] array of coordinates, same as in X
%
%  Outputs
% ---------
% Idx:    [n x 1] indices of X in Y that are closest based on
%                 min(abs(Function(X,Y)))
% IdxVal: [n x 1] difference / distance between X and Y at matching pair
%                 based on min(abs(Function(X,Y)))
%
%  Example usage
% ---------------
% [Idx, IdxVal] = match_min( XYZ1, XYZ2, @pdist2 )
%
%-------------------------------------------------------------------------
% See also MATCH_VECTORS, MATCH2BINS, MATCH2BINS_BEFORE, MATCH2BINS_AFTER
%
%-------------------------------------------------------------------------
% NB: this is essentially like match2bins.m but generalized to other
% (user-defined) functions, like e.g. euclidean distances.
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, July 2019
%-------------------------------------------------------------------------

Idx = nan(size(X,1),1);
IdxVal = nan(size(X,1),1);
for i = 1:length(X)
    [IdxValTemp,IdxTemp] = min(abs(feval(Function,X(i,:),Y)));
    Idx(i) = IdxTemp;
    IdxVal(i) = IdxValTemp;
end

end

